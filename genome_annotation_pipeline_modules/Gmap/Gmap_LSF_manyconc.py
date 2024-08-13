# by Guna
"""
There are two ways to use this file:

The module/class way:
Read __main__ section for example.
Also almost all the properties can be set either via __init__, or via setters.
The demo in __main__ uses __init__

The standalone way:
use the flag -h for help.
Example:
python /opt/common/compbio/lib/python2.7/site-packages/phi/Analyses/Gmap/Gmap_LSF_manyconc.py \
-i /anno/gsap/benchmark/ATH_EST_sequences_20101108.fas \
-d /anno/gsap/benchmark/TAIR10_chr_all.fa \
-o /gsap/dev/data/guna_script_test/output/results.gff3 -t EST \
-D /ngsprod/gsap/GMAPDB -db TAIR10_chr_all \
-l /gsap/dev/data/guna_script_test/output/log.txt -P gtdi-gsap-gmap \
-j 50 -n 50 -J GmapTest -ap='-n 20 -K 50000 --min-identity=0.95 \
--min-trimmed-coverage=0.90' -lp='-R "rusage[mem=2500,scr=100]"' -s -an GMAP
"""

# import some modules and constants
import os
import sys
import time
import re
import argparse
import textwrap
from subprocess import PIPE, Popen
# other Pioneer developed modules and constants
import phi.Logger
import phi.Utils
import phi.Analyses.LSF.Analysis
from phi.Analyses.LSF.Constants import LSF_DEFAULT_QUEUE, \
    LSF_MAX_CONCURRENT_JOB_TOTAL, LSF_DEFAULT_OUTPUT_DIR
from phi.Analyses.Constants import EMAIL_SENDER, EXT_STDOUT, EXT_STDERR, BIN_SET_PIPEFAIL
from phi.LSF.Constants import DEFAULT_MAX_DISK_CHECK_FREQUENCY, \
    DEFAULT_MIN_FREE_SPACE_TO_WARN, DEFAULT_MIN_FREE_SPACE_TO_PAUSE
import phi.LSF.Managers.ManyConcurrent
import phi.DiskSpaceWarning
from .Gmap_Constants import GMAP_ANALYSIS_NAME, \
    GMAP_ANALYSIS_PARAMETERS, LSF_DIR_GMAPDB, \
    GMAP_INDEX_FASTA, GMAP_INDEX_FAS, GMAP_INDEX_FA, \
    GMAP_INDEX_DBFILES, GMAP_TRANSCRIPT_TYPE, GMAP_GFF3_FIELD2, \
    GMAP_OUTPUT_FORMAT_GFF3, LSF_GMAP_PROJECT_NAME, LSF_MAX_JOB_TOTAL, \
    LSF_DELAY_TIME_GMAP, LSF_RESOURCES_GMAP, GMAP_INDEX_GFF3, \
    LSF_BIN_GMAP_V1, LSF_BIN_GMAP_BUILD_V1, \
    LSF_BIN_GMAP_V2, LSF_BIN_GMAP_BUILD_V2, LOCAL_BIN_GMAP_V2
from phi.Parse import cleanFastaNSplit

# Setting Global variables for GMAP program paths - default V1
LSF_BIN_GMAP = LSF_BIN_GMAP_V1
LSF_BIN_GMAP_BUILD = LSF_BIN_GMAP_BUILD_V1

class Gmap_LSF_manyconc(phi.Analyses.LSF.Analysis.Analysis):
    """ Basically inheriting/extending two classes:
1) LSF.Analysis through super class and
2) LSF.managers.ManyConcurrent by creating an object
#1 Not using so many features at this point, but #2 is the main LSF module to submit jobs"""

    def __init__(
            self,
            name=GMAP_ANALYSIS_NAME,
            analysisParameters=GMAP_ANALYSIS_PARAMETERS,
            dbName='',
            dbDestDir=LSF_DIR_GMAPDB, dbFile='',
            inputFile='',
            transcriptType=GMAP_TRANSCRIPT_TYPE[-1],
            outputFile='',
            outputFormat=GMAP_OUTPUT_FORMAT_GFF3[0],
            outputDir='',
            outputDiskName='',
            interimOutputDir='',
            logFh=sys.stderr,
            needEmail=True,
            emails=[phi.Utils.getInteralEmailAddress()],
            verboseLevel=1,
            queue=LSF_DEFAULT_QUEUE,
            projectName=LSF_GMAP_PROJECT_NAME,
            jobName=GMAP_ANALYSIS_NAME,
            lsfParameters='',
            needConcatenateOutput=True,
            jobTotal=LSF_MAX_JOB_TOTAL,
            concurrentJobTotal=LSF_MAX_CONCURRENT_JOB_TOTAL,
            checkExistingStdoutFiles=True,
            requeueable=False
        ):
        #######################################################
        # must declear/define it before setting the super class.
        # Else cannot set lsf manager properties when setting this object's properties.
        self.lsfManager = phi.LSF.Managers.ManyConcurrent.ManyConcurrent()
        # set it in initializer so that user can fine-control it if they want to.
        #######################################################

        # The Analysis super class takes a 2D array for both input and output files
        # (in case some analyses need different kind of input/output
        # files and then they are also split)
        # Currently, inputFilesArray=[[]], outputFilesArray=[[]] are sent empty,
        # but set later after splitting.
        # no filteringParameters='', compressType='', concatenatedOutputFiles=[]
        # no status=None, startTime=None, endTime=None
        phi.Analyses.LSF.Analysis.Analysis.__init__(
            self, name, LSF_BIN_GMAP, analysisParameters, '', dbFile, [[]], [[]],
            outputFormat.lower(), outputDir, outputDiskName, '', needEmail, emails,
            None, None, None, logFh, verboseLevel, queue, projectName, jobName,
            lsfParameters, [], concurrentJobTotal, checkExistingStdoutFiles)

        # gmap specific properties - child class ONLY.
        # Parent's properties are set by above init() call
        # input sequence transcript type - used for GFF column2 (source)
        self.transcriptType = transcriptType
        self.dbName = dbName  # The name of the GMAP DB
        self.dbDestDir = dbDestDir or LSF_DIR_GMAPDB # Common dir for all GMAP databases
        self.inputFile = inputFile  # original input file
        self.inputFiles = [] # smaller child files.
        # Super class sets this too based on inputFile Array size. But here,
        # splitting hasn't happened yet. So, value set by user or default in GMAP_Constants file.
        self.jobTotal = jobTotal or LSF_MAX_JOB_TOTAL
        self.outputFile = outputFile  # the main outputFile
        self.outputFiles = []  # the child output files, can be transient.
        # If interim dir not given by user,
        # it is created at LSF_DEFAULT_OUTPUT_DIR = '/analysis-biocomp02/lsfManager'
        childInputFileDir = LSF_DEFAULT_OUTPUT_DIR + '/' + \
            phi.Utils.getUniqueStringFromHostnameTimePid()
        self.interimOutputDir = interimOutputDir or childInputFileDir
        self.needConcatenateOutput = needConcatenateOutput # this module, it is single one.
        self.lsfDelayTime = LSF_DELAY_TIME_GMAP # for testing, use 30.
        self.jobNames = []    # Will be set after splitting files
        self.commands = []
        self.gff3files = []    # only useful for filtering the raw output to the gff3 files.
        self.toDeleteFiles = []    # Add any temporary files to this array, that are SURE to go.
        self.logger = phi.Logger.Logger(name, logFh)

        # for LSF manager
        self.lsfManager.logFh = logFh
        self.lsfManager.queue = queue
        self.lsfManager.projectName = projectName
        self.lsfManager.jobNames = self.jobNames # Will be set after splitting files
        self.lsfManager.commands = []
        self.lsfManager.requeueable = requeueable
        # Super class value should be over-ridden above
        # by user value or default. Updated after splitting.
        self.lsfManager.jobTotal = self.jobTotal
        self.lsfManager.checkExistingStdoutFiles = checkExistingStdoutFiles
        self.lsfManager.bsubParameters = lsfParameters or LSF_RESOURCES_GMAP
        self.lsfManager.verboseLevel = verboseLevel
        self.lsfManager.concurrentJobTotal = concurrentJobTotal

    # overwrite the default behavior of some setters.
    def __setattr__(self, n, v):
        # aType = type(v).__name__
        #if n == 'dbDestDir':
        #    v = v.lower()

        # if user use attribute way, we need make sure that lsfManager's properties are set too.
        if n == 'logFh' or n == 'queue' or n == 'projectName' or n == 'jobName' or \
            n == 'requeueable' or n == 'jobTotal' or n == 'checkExistingStdoutFiles' or \
            n == 'bsubParameters' or n == 'verboseLevel' or n == 'minFreeSpaceToWarn' or \
            n == 'minFreeSpaceToPause':
            self.lsfManager.__setattr__(n, v)

        phi.Analyses.LSF.Analysis.Analysis.__setattr__(self, n, v)


    # check inputs to make sure they are all ok before trying to run.
    def checkInputs(self):
        """check Inputs"""

        # check on input file - required
        if not self.inputFile or not os.path.exists(self.inputFile):
            self.logger.error("Input fasta file not set or does not exist: %s" % self.inputFile)
            raise Exception, "Input fasta file not set or does not exist: %s" % self.inputFile
        self.inputFile = os.path.abspath(self.inputFile)

        # check on db file - required
        if not self.dbFile or not os.path.exists(self.dbFile):
            self.logger.error("DB fasta file not set or does not exist: %s" % self.dbFile)
            raise Exception, "DB fasta file not set or does not exist: %s" % self.dbFile
        self.dbFile = os.path.abspath(self.dbFile)

        # check on outputDir
        if not os.path.exists(self.outputDir):
            self.logger.error("Output dir does not exist: %s" % self.outputDir)
            raise Exception, "Output dir does not exist: %s" % self.outputDir

        if not os.access(self.outputDir, os.W_OK) or not os.access(self.outputDir, os.X_OK):
            self.logger.error("Output dir %s is either not writable or executable. \
I need both to be able to write to it" % self.outputDir)
            raise Exception, "Output dir %s is either not writable or executable. \
I need both to be able to write to it" % self.outputDir

        # check on outputFile
        if (not self.outputFile) and self.needConcatenateOutput:
            self.outputFile = self.outputDir + '/' + 'results' + GMAP_INDEX_GFF3

        # reset the interimDir to outputDir if no concatenation or non-gff3 output.
        if (not self.needConcatenateOutput) or (self.outputFormat not in GMAP_OUTPUT_FORMAT_GFF3):
            self.interimOutputDir = self.outputDir

        # Interim output dir - parsed & split files are created here - On cluster temp space
        if not os.path.exists(self.interimOutputDir):
            os.mkdir(self.interimOutputDir)

        # -D option for GMAP
        dbDestDir = self.dbDestDir
        if not os.path.exists(dbDestDir):
            self.logger.error("GMAP DB dest dir does not exist: %s" % dbDestDir)
            raise Exception, "GMAP DB dest dir does not exist: %s" % dbDestDir
        self.dbDestDir = os.path.abspath(dbDestDir)

        # -d option for GMAP
        if not self.dbName:
            (dbFileDir, dbFileWithoutDir) = os.path.split(self.dbFile)
            self.dbName = re.sub(GMAP_INDEX_FASTA, '', dbFileWithoutDir, flags=re.IGNORECASE)
            self.dbName = re.sub(GMAP_INDEX_FAS, '', self.dbName, flags=re.IGNORECASE)
            self.dbName = re.sub(GMAP_INDEX_FA, '', self.dbName, flags=re.IGNORECASE)

        # -D /ngsprod/gsap/GMAPDB -d TAIR10_chr_all
        # DB Index files in /ngsprod/gsap/GMAPDB/TAIR10_chr_all:
        # TAIR10_chr_all.chromosome.iit, TAIR10_chr_all.contig.iit ...
        # Checking if all of them exists; Else create the GMAP DB
        dbIndexFilesDir = self.dbDestDir + '/' + self.dbName
        dbIndexFiles = [dbIndexFilesDir + '/' + self.dbName + indx for indx in GMAP_INDEX_DBFILES]

        if all([os.path.exists(dbIndexFile) for dbIndexFile in dbIndexFiles]):
            self.logger.info("Using existing GMAP DB already built at: %s" % dbIndexFilesDir)
        else:
            self.logger.info("GMAP DB does not exist or missing index files. \
Need to create first at: %s" % dbIndexFilesDir)
            self.logger.info("Parsing DB fasta file: %s..." % self.dbFile)
            # Parse DB fasta file to clean; No splitting -
            # returns an array still. Empty if invalid sequences
	    seqNameFile = self.dbDestDir + '/' + self.dbName + '_origSeqNames'
            newDbFiles = cleanFastaNSplit(
                self.dbFile, 1, seqNameFile, 0, self.interimOutputDir, self.dbName, self.logger)
            if len(newDbFiles) == 0:
                self.logger.error("DB fasta file has no valid sequences: %s" % self.dbFile)
                raise Exception, "DB fasta file has no valid sequences: %s" % self.dbFile
            # Run GMAP build
            # gmap_build -d TAIR10_chr_all -D /ngsprod/gsap/GMAPDB
            # /anno/gsap/benchmark/TAIR10_chr_all.fa
            command = ' '.join([
                LSF_BIN_GMAP_BUILD, '-d', self.dbName, '-D', self.dbDestDir, newDbFiles[0]])
            self.logger.info("Running GMAP DB build command: %s" % command)
            returnCode = os.system(command)
            if returnCode != 0:
                self.logger.error("DB build Failed. command %s did not return 0: %d" % \
                    (command, returnCode))
                raise Exception, "DB build Failed. command %s did not return 0: %d" % \
                    (command, returnCode)
            self.logger.info("GMAP DB build done: %s" % dbIndexFilesDir)
            self.dbFile = newDbFiles[0]
            self.toDeleteFiles.append(newDbFiles[0])


        # Parsing & Splitting the input file
        if self.verboseLevel > 0:
            self.logger.info("Trying to split input seq file %s up to %d \
files at %s" % (self.inputFile, self.jobTotal, self.interimOutputDir))
        self.inputFiles = cleanFastaNSplit(
            self.inputFile, self.jobTotal, None, 0, self.interimOutputDir, None, self.logger)


        if len(self.inputFiles) == 0:
            self.logger.error("Input fasta file has no valid sequences: %s" % self.inputFile)
            raise Exception, "Input fasta file has no valid sequences: %s" % self.inputFile
        for f in self.inputFiles:
            if not os.path.exists(f):
                self.logger.error("Input file does not exist: %s" % f)
                raise Exception, "Input file does not exist: %s" % f
                # need to test readable or not, but there is no easy way.
        # adding to super class attribute, so it could be used for delete functions later
        self.inputFilesArray.append(self.inputFiles)

        # In case fewer files are created, than max LSF job total
        self.jobTotal = len(self.inputFiles)
        self.lsfManager.jobTotal = self.jobTotal
        if self.verboseLevel > 0:
            self.logger.info("Input seq file is actually split into only \
%d files at %s" % (self.jobTotal, self.interimOutputDir))

        #self.jobNames = []
        for i in range(0, self.jobTotal):
            fileNumber = str(i + 1)
            self.jobNames.append(self.jobName + '_' + fileNumber)
        self.lsfManager.jobNames = self.jobNames

        #check if valid output format
        if self.outputFormat not in GMAP_OUTPUT_FORMAT_GFF3:
            self.logger.error("Invalid output format entered for GMAP: %s. \
Should be one of the gff formats: %s" % (self.outputFormat, GMAP_OUTPUT_FORMAT_GFF3))
            raise Exception, "Invalid output format entered for GMAP: %s. \
Should be one of the gff formats: %s" % (self.outputFormat, GMAP_OUTPUT_FORMAT_GFF3)

        return


    def run(self):
        """run the analysis"""

        self.checkInputs()

        # give SONAS some time to let the newly created file to appear on the cluster nodes.
        if self.verboseLevel > 0:
            self.logger.info("give SONAS %d seconds to let the newly created \
files to appear on the cluster nodes" % self.lsfDelayTime)
        time.sleep(self.lsfDelayTime)

        # check space:
        # causing python subprocess os.fork() memory issues; too much trouble for a minor benign check
        #phi.DiskSpaceWarning.checkDiskSpace(
        #    self.outputDiskName, self.minFreeSpaceToWarn,
        #    self.minFreeSpaceToPause, self.name, EMAIL_SENDER, self.emails,
        #    DEFAULT_MAX_DISK_CHECK_FREQUENCY, self.logFh)

        # create lsf manager obj.
        if self.verboseLevel > 0:
            self.logger.info("creating lsf manager...")

        lsfManager = self.lsfManager
        # make command stdout stderr etc.
        countInput = 0
        for inputFile in self.inputFiles:
            # Only one GFF3 output file.
            # Allowed output formats: GMAP_OUTPUT_FORMAT_GFF3 =
            # ['gff3_gene', 'gff3_match_cdna', 'gff3_match_est']
            (inputDir, inputFileWithoutDir) = os.path.split(inputFile)
            outputFile = inputDir + '/' + inputFileWithoutDir.replace(GMAP_INDEX_FA, '') + \
                self.outputFormat.replace('gff3', '') + GMAP_INDEX_GFF3
            self.outputFiles.append(outputFile)

            # merging only gff3 files. In this case, gff3 file = output file
            gff3File = outputFile
            self.gff3files.append(gff3File)

            countInput += 1
            stdoutFile = self.outputDir + '/' + str(countInput) + EXT_STDOUT
            stderrFile = self.outputDir + '/' + str(countInput) + EXT_STDERR
            command = ''

            # gmap -d TAIR10_chr_all -D /ngsprod/gsap/GMAPDB -f gff3_gene
            # -n 20 -K 50000 --min-identity=0.95 --min-trimmed-coverage=0.90
            # /anno/gsap/benchmark/ATH_EST_sequences_20101108.fas > AT_EST_TAIR10.gff3
            command = ' '.join([
                BIN_SET_PIPEFAIL,
                'cd ' + self.interimOutputDir,
                '&&',
                LSF_BIN_GMAP,
                '-d',
                self.dbName,
                '-D',
                self.dbDestDir,
                '-f',
                self.outputFormat,
                self.analysisParameters,
                inputFile,
                '>',
                gff3File])

            self.logger.info(command)
            self.commands.append(command)
            self.stdoutFiles.append(stdoutFile)
            self.stderrFiles.append(stderrFile)

        self.outputFilesArray.append(self.outputFiles)

        # sub LSF jobs.
        lsfManager.submit(self.commands, lsfManager.jobNames, self.stdoutFiles, self.stderrFiles)
        # run LSF jobs.
        lsfManager.checkManyJobsUntilAllDone()
        if self.verboseLevel > 0:
            self.logger.info("give SONAS %d seconds to let the newly created \
output files to appear on the submission host." % self.lsfDelayTime)
        time.sleep(self.lsfDelayTime)

        # print out stats.
        self.logsObj = lsfManager.setLogsObj()
        self.logger.info("Cluster stats:")
        self.logger.info("CPU total: %d sec." % self.logsObj.sumCpuTime)
        self.logger.info("CPU max: %d sec." % self.logsObj.maxCpuTime)
        self.logger.info("CPU min: %d sec." % self.logsObj.minCpuTime)
        self.logger.info("CPU ave: %s sec." % self.logsObj.aveCpuTimeStr)
        self.logger.info("Memory max: %d MB" % self.logsObj.maxMemory)
        self.logger.info("Memory min: %d MB" % self.logsObj.minMemory)
        self.logger.info("Memory ave: %s MB" % self.logsObj.aveMemoryStr)

        # concatenation, only for gff3. else I cannot concatenate the raw binary files
        if self.needConcatenateOutput and (self.outputFormat in GMAP_OUTPUT_FORMAT_GFF3):
            # check space as the output size will be doubled temporarily.
            phi.DiskSpaceWarning.checkDiskSpace(
                self.outputDiskName, self.minFreeSpaceToWarn,
                self.minFreeSpaceToPause, self.name, EMAIL_SENDER,
                self.emails, DEFAULT_MAX_DISK_CHECK_FREQUENCY, self.logFh)

            if self.verboseLevel > 0:
                self.logger.info("Trying to concatenate the outputs into \
one file %s at %s" % (self.outputFile, phi.Utils.getTimeStampString()))

            if self.checkExistingStdoutFiles:
                # first test if the concatenated files have been
                # created before and child files are deleted.
                missingStdoutFileTotal = 0
                for f in self.stdoutFiles:
                    if not os.path.exists(f):
                        missingStdoutFileTotal += 1

                if missingStdoutFileTotal == 0 and os.path.exists(self.outputFile):
                    if self.verboseLevel > 0:
                        self.logger.info("Looks like concatenation was \
done before. All stdout files look ok. The output file %s is present. \
No rerun at %s" % (self.outputFile, phi.Utils.getTimeStampString()))
                    self.endTime = time.localtime()
                    return


            fhOut = open(self.outputFile, 'w')

            # not simple
            for f in self.gff3files:
                # if no hit for gmap, no output at all.
                if not os.path.exists(f):
                    continue

                fhIn = open(f)
                for line in fhIn:
                    if line[0] == '#':
                        # comment lines: just copy it.
                        fhOut.write(line)
                    else:
                        # feature lines: convert field 2 to GMAP_transcriptType
                        row = line.split('\t')
                        #refSeqId = row[0]
                        #attrs = row[8]
                        #attrs = attrs.replace('ID=', 'ID=' + refSeqId + '_')
                        #attrs = attrs.replace('Parent=', 'Parent=' + refSeqId + '_')
                        source = GMAP_GFF3_FIELD2 + '_' + self.transcriptType
                        newLine = '\t'.join([
                            row[0], source, row[2], row[3], row[4], row[5], row[6], row[7], row[8]
                        ])
                        fhOut.write(newLine)
                fhIn.close()
            fhOut.close()

        self.endTime = time.ctime()
        self.logger.info("Analysis All well done at %s" % self.endTime)

        self.postProcess()

    def postProcess(self):
        """post process"""

        if self.needEmail:
            logsObj = self.logsObj
            subject = 'Analysis %s with LSF project %s is well done at %s' % \
                (self.name, self.lsfManager.projectName, self.endTime)
            files = []
            logFile = self.lsfManager.logFh.name
            if self.needConcatenateOutput and (self.outputFormat in GMAP_OUTPUT_FORMAT_GFF3):
                files = [self.outputFile, logFile] # + self.lsfManager.stdoutFiles
            else:
                files = [logFile] # + self.lsfManager.stdoutFiles

            # somehow, outlook does not honor \n for the first a few lines.
            body = "Your command (quotes might have been lost): %s\n" % self.commandLine
            body = body + "Actual total jobs: %d\n" % self.jobTotal
            body = body + "Concurrent jobs: %d\n" % self.concurrentJobTotal
            body = body + "Cluster stats:\n"
            body = body + "CPU total: %d sec.\n" % logsObj.sumCpuTime
            body = body + "CPU max: %d sec.\n"   % logsObj.maxCpuTime
            body = body + "CPU min: %d sec.\n"   % logsObj.minCpuTime
            body = body + "CPU ave: %s sec.\n"   % logsObj.aveCpuTimeStr
            body = body + "Memory max: %d MB\n"  % logsObj.maxMemory
            body = body + "Memory min: %d MB\n"  % logsObj.minMemory
            body = body + "Memory ave: %s MB\n"  % logsObj.aveMemoryStr
            body = body + 'For more details, please check the files under %s\n' % \
                self.outputDir
            body = body + 'Please also read the following result file(s) and \
log file:\n' + '\n'.join(files)

            for email in self.emails:
                phi.Utils.sendEmail(EMAIL_SENDER, email, subject, body)

        # Delete any temporary files (Ex: Newly created DB file)
        # - this function expects an array of files
        self.deleteFiles(self.toDeleteFiles)
        # delete child input files
        self.deleteInputFiles()
        # now delete the child output files so as to save the space.
        if self.needConcatenateOutput and (self.outputFormat in GMAP_OUTPUT_FORMAT_GFF3):
            self.deleteOutputFiles()

        # This would be redundant here
        if self.needConcatenateOutput:
            self.deleteFiles(self.gff3files)

        # Deleting temp files on Wengang's request
        self.deleteStderrFiles()
        self.deleteStdoutFiles()
        phi.Utils.safermtree(self.interimOutputDir)
        self.logger.info("Deleting temp data well done at %s" % time.ctime())

# Help Menu formatting
# Trivial new class to maintain 2 different formatting in the command-line help menu below
class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

# Helper function to return stdout from a command
def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True)
    #process.communicate() returns an array [stdout, stderr]
    return process.communicate()[0]

# Custom action (inspired by _VersionAction) to output analysis program's version info
# There could be multiple versions for the analysis program. ex: gmap
class CustomAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 program=None,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help="show analysis program's version and exit"):
        super(CustomAction, self).__init__(
		option_strings=option_strings,
		dest=dest,
		default=default,
		nargs=0,
		help=help)
        self.program = program
    def __call__(self, parser, namespace, values, option_string=None):
        program = self.program
        if (isinstance(program, list)):
            for i, prog in enumerate(program):
                print("Version %d - analysis program invoked: %s" %(i+1, prog))
                print(cmdline(prog + ' --version'))
        else:
	    print("Analysis program invoked: %s" %(program))
	    print(cmdline(program + ' --version'))
        #message = cmdline(program + ' --version')
        #formatter = parser._get_formatter()
        #formatter.add_text(message)
        #parser._print_message(formatter.format_help(), argparse._sys.stdout)
        parser.exit()


def run():
    """run analysis"""

    parser = argparse.ArgumentParser(
        description="Command line tool to run GMAP using HPC cluster",
        epilog=textwrap.dedent('''\
Not all combinations of parameters are supported. User is responsible for NOT passing conflicting parameters into the script.
                                  
Example:
python %(prog)s -i /anno/gsap/benchmark/ATH_EST_sequences_20101108.fas -d 
/anno/gsap/benchmark/TAIR10_chr_all.fa 
-o /gsap/dev/data/guna_script_test/output/results.gff3 -t EST 
-D /ngsprod/gsap/GMAPDB -db TAIR10_chr_all
-l /gsap/dev/data/guna_script_test/output/log.txt 
-P gtdi-gsap-gmap -j 50 -n 50 -J GmapTest -ap='-n 20 -K 50000
--min-identity=0.95 --min-trimmed-coverage=0.90'
-lp='-R "rusage[mem=2500,scr=100]"' -s -an GMAP
'''),
        formatter_class=CustomFormatter)

    # Mandatory arguments
    parser.add_argument(
        '-i', dest='inputFile', required=True,
        help="**REQUIRED** input fasta file, DNA. It does not have \
to be on a cluster-accessible file system.")
    parser.add_argument(
        '-d', dest='dbFile', required=True,
        help="**REQUIRED** database file, DNA. If GMAP database doesn't exist, \
 the program will try to build (-D, -db options). The filename part should be \
 xxx.fa. xxx part should be alphanumeric.")
    parser.add_argument(
        '-o', dest='outputFile', required=True,
        help="**REQUIRED** output dir (unconcatenated output) or output file \
(concatenated gff3 output file). Full path. Must be cluster accessible.")

    # Optional arguments
    parser.add_argument(
        '-t', dest='transcriptType', default=GMAP_TRANSCRIPT_TYPE[-1],
        help="input sequence transcript type. Ex: EST, cDNA, Assembly, \
etc. Shows in GFF3 annotation (source column).")
    parser.add_argument(
        '-D', dest='dbDestDir', default=LSF_DIR_GMAPDB,
        help="destination dir for GMAP database. If provided, dir should \
exist. GMAP DB for the dbFile will be built here if missing index files.")
    parser.add_argument(
        '-db', dest='dbName',
        help="GMAP database name (sub-directory within dbDestDir). \
Will be prefix of dbFile if not given.")
    parser.add_argument(
        '-an', dest='analysisName', default=GMAP_ANALYSIS_NAME,
        help="analysis name. Used in email report.")
    parser.add_argument(
        '-of', dest='outputFormat', choices=GMAP_OUTPUT_FORMAT_GFF3,
        default=GMAP_OUTPUT_FORMAT_GFF3[0],
        help="output format. currently only GFF3 formats possible.")
    parser.add_argument(
        '-ap', dest='analysisParameters', default=GMAP_ANALYSIS_PARAMETERS,
        help="GMAP specific analysis parameters. Provide as one string within \
quotes with equal-to sign. Ex: -ap='-n 20 -K 50000 --min-identity=0.95 \
--min-trimmed-coverage=0.90'. Does not conflict with other program-specific options.")
    parser.add_argument(
        '-lp', dest='lsfParameters', default=LSF_RESOURCES_GMAP,
        help='''lsf parameters. Provide as one string within quotes \
with equal-to sign. Ex: -lp='-R "rusage[mem=2500,scr=100]"'. \
For special options only such as resources, not for queue, \
project name, job name etc. Do not use it if you are not sure.''')
    parser.add_argument(
        '-ld', dest='lsfDelayTime', type=float, default=LSF_DELAY_TIME_GMAP,
        help="LSF delay time in seconds. Wait this long before \
LSF job submission and after LSF job completion.")
    parser.add_argument(
        '-q', dest='queue', default=LSF_DEFAULT_QUEUE, help="lsf queue name")
    parser.add_argument(
        '-P', dest='projectName', default=LSF_GMAP_PROJECT_NAME,
        help="lsf project name")
    parser.add_argument(
        '-J', dest='jobName', default=GMAP_ANALYSIS_NAME, help="lsf job name")
    parser.add_argument(
        '-j', dest='maxLsfJob', type=int, default=LSF_MAX_JOB_TOTAL,
        help="max LSF jobs.")
    parser.add_argument(
        '-n', dest='maxConcurrentLsfJob', type=int,
        default=LSF_MAX_CONCURRENT_JOB_TOTAL,
        help="max concurrent LSF jobs.")
    parser.add_argument(
        '-l', dest='logFile', help="log file (not for LSF jobs), full path. \
Default: stderr. It does not have to be on a cluster-accessible file system.")
    parser.add_argument(
        '-s', '--needEmail', action='store_true', help="Flag set to send email \
upon being successfully done. Default: no emails. Email will always be sent \
out in cases of errors.")
    parser.add_argument(
        '-e', dest='emails', help="email addresses. Delimited by space \
if multiple. eg: -e='guna.gurazada@pioneer.com another@pioneer.com'. \
Default: submitter")
    parser.add_argument(
        '-z', dest='bsubInterval', type=float, default=1.0,
        help="wait time in seconds, between bsub commands. \
Can be a float, but must be less than 10.")
    parser.add_argument(
        '-r', '--requeueable', action='store_true', help="Flag set to requeue \
jobs. Default: not. If there is still a problem after 2 tries, the status will \
become PSUSP, you need to manually kill both the process and LSF jobs. However, \
there are other reasons that a status become PSUSP. Check LSF stdout and stderr to make sure.")
    parser.add_argument(
        '-m', dest='minFreeSpaceToWarn', type=int,
        default=DEFAULT_MIN_FREE_SPACE_TO_WARN/1000000,
        help="min free disk space to warn. Integer in MB.")
    parser.add_argument(
        '-M', dest='minFreeSpaceToPause', type=int,
        default=DEFAULT_MIN_FREE_SPACE_TO_PAUSE/1000000,
        help="min free disk space to pause. Integer in MB.")
    parser.add_argument(
        '-S', '--checkExistingStdoutFiles', action='store_true',
        help="check previous LSF stdout files. Used in case of \
resubmit to avoid the redundant LSF jobs.")
    parser.add_argument(
        '-v', dest='verboseLevel', type=int, choices=[0, 1, 2], default=1,
        help="verbose level. 0: quiet; 1: important messages; 2: all messages. \
Errors will always be spit out regardless the setting.")
    parser.add_argument(
        '-V', '--useVersion', dest='gmapVersion', type=int, choices=[1, 2], default=1,
        help="GMAP version to run. Use --version flag to see which GMAP versions are installed \
and set this option accordingly. Make sure to use appropriate -D dbDestDir and -ap analysisParameters \
options for the version selected.")
    parser.add_argument(
        '--version', action=CustomAction, program=[LSF_BIN_GMAP_V1, LOCAL_BIN_GMAP_V2],
        help="show analysis program's version and exit")

    
    args = parser.parse_args()

    # new parameters
    maxBsubInterval = 10.0
    childInputFileDir = LSF_DEFAULT_OUTPUT_DIR + '/' + \
        phi.Utils.getUniqueStringFromHostnameTimePid()  # Interim temp dir

    # open log file first
    logFh = None
    if not args.logFile:
        logFh = sys.stderr
    else:
        logFh = open(args.logFile, 'a', 0) # no buffering.

    # parsing options
    inputFile = args.inputFile
    dbFile = args.dbFile
    outputFile = args.outputFile
    if outputFile[0] != '/':
        sys.stderr.write("output file %s in not in full path\n" % outputFile)
        logFh.write("output file %s in not in full path\n" % outputFile)
        exit(1)

    if os.path.isdir(outputFile):
        outputDir = outputFile
        outputFile = ''
    else:
        outputDir = os.path.dirname(outputFile)

    transcriptType = args.transcriptType
    dbDestDir = args.dbDestDir
    dbName = args.dbName
    analysisName = args.analysisName
    outputFormat = args.outputFormat
    analysisParameters = args.analysisParameters
    lsfParameters = args.lsfParameters
    lsfDelayTime = args.lsfDelayTime
    queue = args.queue
    projectName = args.projectName
    jobName = args.jobName
    needConcatenateOutput = True  # Removed this as user option & always default to True
    maxLsfJob = args.maxLsfJob
    maxConcurrentLsfJob = args.maxConcurrentLsfJob
    needEmail = args.needEmail
    if not args.emails:
        emails = [phi.Utils.getInteralEmailAddress()]
    else:
        emails = args.emails.split(' ')
    bsubInterval = args.bsubInterval
    requeueable = args.requeueable
    minFreeSpaceToWarn = args.minFreeSpaceToWarn * 1000000  # input is in mega
    minFreeSpaceToPause = args.minFreeSpaceToPause * 1000000 # input is in mega
    checkExistingStdoutFiles = args.checkExistingStdoutFiles
    verboseLevel = args.verboseLevel
    gmapVersion = args.gmapVersion
    
    # change GMAP related paths based on the GMAP version provided
    # defaults are set for version 1 - so only change if version 2 is given 
    if gmapVersion == 2:
	global LSF_BIN_GMAP, LSF_BIN_GMAP_BUILD
	LSF_BIN_GMAP = LSF_BIN_GMAP_V2
	LSF_BIN_GMAP_BUILD = LSF_BIN_GMAP_BUILD_V2

    # check on a few options before submitting in the try block
    if len(projectName) < 8:
        sys.stderr.write("project name: %s is too short or in wrong format.\n" % projectName)
        logFh.write("project name: %s is too short or in wrong format.\n" % projectName)
        exit(1)
    if bsubInterval > maxBsubInterval:
        sys.stderr.write("bsubInterval %f is greater than max %f.\n" % \
            (bsubInterval, maxBsubInterval))
        logFh.write("bsubInterval %f is greater than max %f.\n" % \
            (bsubInterval, maxBsubInterval))
        exit(1)
    if verboseLevel > 0:
        logFh.write("The command is as below (quotes might have been removed \
when restoring the command):\n%s\n" % ' '.join(sys.argv))

    # now comes the time-consuming part. Use try to catch the errors and report errors.
    # there are too many types of errors and impossible to check individually
    # (eg: network, space, hard drive failure). Just trace all
    try:
        # now create the object. If you are writing your own script,
        # you need to import it and use the full name:
        # phi.Analyses.Gmap.Gmap_LSF_manyconc.Gmap_LSF_manyconc()
        analysisObj = Gmap_LSF_manyconc(
            analysisName, analysisParameters, dbName, dbDestDir, dbFile, inputFile,
            transcriptType, outputFile, outputFormat, outputDir, outputDir,
            childInputFileDir, logFh, needEmail, emails, verboseLevel,
            queue, projectName, jobName, lsfParameters, needConcatenateOutput,
            maxLsfJob, maxConcurrentLsfJob, checkExistingStdoutFiles, requeueable)
        analysisObj.minFreeSpaceToWarn = minFreeSpaceToWarn
        analysisObj.minFreeSpaceToPause = minFreeSpaceToPause
        analysisObj.lsfDelayTime = lsfDelayTime
        analysisObj.lsfManager.bsubInterval = bsubInterval

        analysisObj.commandLine = ' '.join(sys.argv)
        analysisObj.run()

    except:
        import traceback
        stackTrace = traceback.format_exc() + \
            "\nSome jobs might be still running and you need manually bkill them\n"
        # log the error
        logFh.write(stackTrace)

        # email the users. For testing purpose, disable emailing.
        for email in emails:
            phi.Utils.sendEmail(EMAIL_SENDER, email, 'Error in ' + analysisName, stackTrace)

        raise
