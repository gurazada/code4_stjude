
"""
There are two ways to use this file:

The module/class way:
Read __main__ section for example.
Also almost all the properties can be set either via __init__, or via setters.
The demo in __main__ uses __init__

The standalone way:
use the flag -h for help.
Example:
repeatmasker_lsf.py -i guna-testing/RepeatMasker/input/ATH_EST_sample_clean_9.fa \
-lib /ngsprod/gsap/genomes/mipsREdat_9.3p_ALL.fasta \
-dir guna-testing/RepeatMasker/output/ATH_sample -o ATH_EST_sample_clean_9 \
-gff -x -l guna-testing/RepeatMasker/output/ATH_sample/log.txt \
-P gtdi-gsap-repeatmasker -j 50 -n 50 -J RMtest -ap='-nolow' \
-lp='-R "rusage[mem=6000,scr=500]"' -s -an RM
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
    LSF_DEFAULT_OUTPUT_DIR, LSF_MAX_CONCURRENT_JOB_TOTAL
from phi.Analyses.Constants import EMAIL_SENDER, \
    DEFAULT_MAX_DISK_CHECK_FREQUENCY, EXT_STDOUT, EXT_STDERR, BIN_SET_PIPEFAIL
from phi.LSF.Constants import DEFAULT_MIN_FREE_SPACE_TO_WARN, \
     DEFAULT_MIN_FREE_SPACE_TO_PAUSE
import phi.LSF.Managers.ManyConcurrent
import phi.DiskSpaceWarning
from .RM_Constants import RM_ANALYSIS_NAME, RM_ANALYSIS_PARAMETERS, \
    LSF_RM_PROJECT_NAME, LSF_MAX_JOB_TOTAL, LSF_BIN_RM, LSF_MAX_JOB_TOTAL, \
    LSF_DELAY_TIME_RM, LSF_RESOURCES_RM, RM_INDEX_FASTA, RM_INDEX_OUT_MASKED, \
    RM_INDEX_OUT_GFF, RM_INDEX_OUTFILES_UNWANTED
from phi.Parse import cleanFastaNSplit, mergeSeqSlices

class RM_LSF_manyconc(phi.Analyses.LSF.Analysis.Analysis):
    """ Basically inheriting/extending two classes:
        1) LSF.Analysis through super class and
        2) LSF.managers.ManyConcurrent by creating an object
        #1 Not using so many features at this point, but #2 is the main LSF module to submit jobs"""

    def __init__(
            self,
            name=RM_ANALYSIS_NAME,
            analysisParameters=RM_ANALYSIS_PARAMETERS,
            inputFile='',
            species='',
            dbFile='',
            outputFilePrefix='',
            outputDir='',
            chunkSize=0,
            outputGff=False,
            maskWithX=False,
            maskWithSmall=False,
            outputDiskName='',
            interimOutputDir='',
            logFh=sys.stderr,
            needEmail=True,
            emails=[phi.Utils.getInteralEmailAddress()],
            verboseLevel=1,
            queue=LSF_DEFAULT_QUEUE,
            projectName=LSF_RM_PROJECT_NAME,
            jobName=RM_ANALYSIS_NAME,
            lsfParameters='',
            jobTotal=LSF_MAX_JOB_TOTAL,
            concurrentJobTotal=LSF_MAX_CONCURRENT_JOB_TOTAL,
            checkExistingStdoutFiles=True,
            requeueable=False
        ):
        #######################################################
        # must declare/define it before setting the super class.
        # Else cannot set lsf manager properties when setting this object's properties.
        # set it in initializer so that user can fine-control it if they want to.
        self.lsfManager = phi.LSF.Managers.ManyConcurrent.ManyConcurrent()
        #######################################################

        # The Analysis super class takes a 2D array for both input and output files
        # (in case some analyses need different kind of
        # input/output files and then they are also split)
        # Currently, inputFilesArray=[[]], outputFilesArray=[[]]
        # are sent empty, but set later after splitting.
        # no filteringParameters='', outputFormat='', compressType='', concatenatedOutputFiles=[]
        # no status=None, startTime=None, endTime=None
        phi.Analyses.LSF.Analysis.Analysis.__init__(
            self, name, LSF_BIN_RM, analysisParameters, '', dbFile, [[]],
            [[]], '', outputDir, outputDiskName, '', needEmail, emails,
            None, None, None, logFh, verboseLevel, queue, projectName,
            jobName, lsfParameters, [], concurrentJobTotal, checkExistingStdoutFiles)

        # RepeatMasker specific properties - child class ONLY.
        # Parent's properties are set by above init() call
        # dbFile, outputDir are set through above init() call
        # split fasta sequences into smaller chunks. default: 0. Merging also changes based on this.
        self.chunkSize = chunkSize
        # -species option for RepeatMasker if using repeat masker's inbuilt libraries
        self.species = species
        # -gff option for RepeatMasker creates additional gff file
        self.outputGff = outputGff
        # -x option for RepeatMasker masks repeats with Xs instead of Ns
        self.maskWithX = maskWithX
         # -xsmall option for RepeatMasker masks repeats with lowercase
        self.maskWithSmall = maskWithSmall
        # will be stripped from inputFile if not given
        self.outputFilePrefix = outputFilePrefix
        self.inputFile = inputFile  # original input file
        self.inputFiles = [] # smaller child files.
        # Super class sets this too based on inputFile Array size. But here,
        # splitting hasn't happened yet. So, value set by user or default in RM_Constants file.
        self.jobTotal = jobTotal or LSF_MAX_JOB_TOTAL
        self.outputFile = ''    # the main outputFile - masked.fa
        self.outputFiles = []     # the child output files, can be transient.
        # If interim dir not given by user, it is created at
        # LSF_DEFAULT_OUTPUT_DIR = '/analysis-biocomp02/lsfManager'
        childInputFileDir = LSF_DEFAULT_OUTPUT_DIR + '/' + \
            phi.Utils.getUniqueStringFromHostnameTimePid()
        self.interimOutputDir = interimOutputDir or childInputFileDir
        self.lsfDelayTime = LSF_DELAY_TIME_RM # for testing, use 30.
        self.jobNames = []    # Will be set after splitting files
        self.commands = []
        self.gff3file = ''     # output gff3 file only if -gff option is given.
        self.gff3files = []    # child gff3 files only if -gff option is given.
        self.toDeleteFiles = []    # Add any temporary files to this array, that are SURE to go.
        self.logger = phi.Logger.Logger(name, logFh)

        # for LSF manager
        self.lsfManager.logFh = logFh
        self.lsfManager.queue = queue
        self.lsfManager.projectName = projectName
        self.lsfManager.jobNames = self.jobNames # Will be set after splitting files
        self.lsfManager.commands = []
        self.lsfManager.requeueable = requeueable
        # Super class value should be over-ridden above by user value
        # or default. Updated after splitting.
        self.lsfManager.jobTotal = self.jobTotal
        self.lsfManager.checkExistingStdoutFiles = checkExistingStdoutFiles
        self.lsfManager.bsubParameters = lsfParameters or LSF_RESOURCES_RM
        self.lsfManager.verboseLevel = verboseLevel
        self.lsfManager.concurrentJobTotal = concurrentJobTotal

    # overwrite the default behavior of some setters.
    def __setattr__(self, n, v):
        # if user use attribute way, we need make sure that lsfManager's properties are set too.
        if n == 'logFh' or n == 'queue' or n == 'projectName' or \
              n == 'jobName' or n == 'requeueable' or n == 'jobTotal' or \
              n == 'checkExistingStdoutFiles' or n == 'bsubParameters' or \
              n == 'verboseLevel' or n == 'minFreeSpaceToWarn' or \
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

	# check on db file (repeat masker library file) and species option - only one should be set, not both
	if self.dbFile and self.species:
	    self.logger.error("Looks like both -species '%s' and -lib '%s' have been set. \
Only one of them is allowed." % (self.species, self.dbFile))
	    raise Exception, "Looks like both -species '%s' and -lib '%s' have been set. \
Only one of them is allowed." % (self.species, self.dbFile)	

        # check on db file (repeat masker library file) and species - one of them is required
        if not self.species:
            if not self.dbFile or not os.path.exists(self.dbFile):
                self.logger.error("Neither species is given, nor repeat masker library file \
is set or exists: %s" % self.dbFile)
                raise Exception, "Neither species is given, nor repeat masker library file \
is set or exists: %s" % self.dbFile
	    else:
		self.dbFile = os.path.abspath(self.dbFile)

        # check on outputDir
        if not os.path.exists(self.outputDir):
            self.logger.warn("Output dir does not exist. Trying to create dir: %s" % self.outputDir)
            try:
                os.mkdir(self.outputDir)
            except OSError:
                self.logger.error("Output dir %s is either invalid path or \
no write permissions." % self.outputDir)
                raise Exception("Output dir %s is either invalid path or \
no write permissions." % self.outputDir)

        if not os.access(self.outputDir, os.W_OK) or not os.access(self.outputDir, os.X_OK):
            self.logger.error("Output dir %s is either not writable or executable. \
I need both to be able to write to it" % self.outputDir)
            raise Exception, "Output dir %s is either not writable or executable. \
I need both to be able to write to it" % self.outputDir
        self.outputDir = os.path.abspath(self.outputDir)

        # check on output file prefix
        if (not self.outputFilePrefix):
            (inputDir, inputFileWithoutDir) = os.path.split(self.inputFile)
            self.outputFilePrefix = inputFileWithoutDir
            for fasta_extn in RM_INDEX_FASTA:
                self.outputFilePrefix = re.sub(
                    fasta_extn, '', self.outputFilePrefix, flags=re.IGNORECASE)


        # output files - final concatenated masked fasta file and gff file - created in outputDir
        self.outputFile = self.outputDir + '/' + self.outputFilePrefix + '_masked.fa'
        if self.outputGff:
            self.gff3file = self.outputDir + '/' + self.outputFilePrefix + '_repeats.gff'


        # Interim output dir - parsed & split files are created here and
        # their respective interim output files - On cluster temp space
        if not os.path.exists(self.interimOutputDir):
            self.logger.warn("Interim Output dir does not exist. Trying to \
create dir: %s" % self.interimOutputDir)
            try:
                os.mkdir(self.interimOutputDir)
            except OSError:
                self.logger.error("Interim Output dir %s is either invalid \
path or no write permissions." % self.interimOutputDir)
                raise Exception("Interim Output dir %s is either invalid path \
or no write permissions." % self.interimOutputDir)

        if not os.access(self.interimOutputDir, os.W_OK) or \
            not os.access(self.interimOutputDir, os.X_OK):
            self.logger.error("Interim Output dir %s is either not writable \
or executable. I need both to be able to write to it" % self.interimOutputDir)
            raise Exception("Interim Output dir %s is either not writable or \
executable. I need both to be able to write to it" % self.interimOutputDir)
        self.interimOutputDir = os.path.abspath(self.interimOutputDir)

        # analysis parameters: -lib, -species, -gff, -x, -xsmall
        if self.dbFile:
            self.analysisParameters = self.analysisParameters + ' -lib ' + self.dbFile
        if self.species:
            self.analysisParameters = self.analysisParameters + ' -species ' + self.species
        if self.outputGff:
            self.analysisParameters = self.analysisParameters + ' -gff'
        if self.maskWithX:
            self.analysisParameters = self.analysisParameters + ' -x'
        if self.maskWithSmall:
            self.analysisParameters = self.analysisParameters + ' -xsmall'


        # check space - in temp dir
        phi.DiskSpaceWarning.checkDiskSpace(
            self.interimOutputDir, self.minFreeSpaceToWarn,
            self.minFreeSpaceToPause, self.name, EMAIL_SENDER,
            self.emails, DEFAULT_MAX_DISK_CHECK_FREQUENCY, self.logFh)

        # Parsing & Splitting the input file
        if self.verboseLevel > 0:
            self.logger.info("Trying to split input seq file %s up to %d \
files at %s" % (self.inputFile, self.jobTotal, self.interimOutputDir))
        self.inputFiles = cleanFastaNSplit(
            self.inputFile, self.jobTotal, None, self.chunkSize, self.interimOutputDir, None, self.logger)


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


        return


    def run(self):
        """run the analysis"""

        self.checkInputs()

        # give SONAS some time to let the newly created file to appear on the cluster nodes.
        if self.verboseLevel > 0:
            self.logger.info("give SONAS %d seconds to let the newly created \
files to appear on the cluster nodes" % self.lsfDelayTime)
        time.sleep(self.lsfDelayTime)

        # check space: in temp dir as RepeatMasker output goes there first
        phi.DiskSpaceWarning.checkDiskSpace(
            self.interimOutputDir, self.minFreeSpaceToWarn,
            self.minFreeSpaceToPause, self.name, EMAIL_SENDER,
            self.emails, DEFAULT_MAX_DISK_CHECK_FREQUENCY, self.logFh)

        # create lsf manager obj.
        if self.verboseLevel > 0:
            self.logger.info("creating lsf manager...")

        lsfManager = self.lsfManager
        # make command stdout stderr etc.
        countInput = 0
        for inputFile in self.inputFiles:
            # RepeatMasker creates a number of output files in the outputDir:
            # .alert, .cat(.gz), .masked, .ori.out, .out, .out.gff, .tbl
            # (appended to input file name)
            # Only keeping .masked as output file & .out.gff as gff file if -gff param is given
            # Just noting the filenames here for merging/deleting later -
            # RepeatMasker only needs the output dir
            # for repeatmasker, Output dir = Interim output dir
            (inputDir, inputFileWithoutDir) = os.path.split(inputFile)
            outputFile = inputDir + '/' + inputFileWithoutDir + RM_INDEX_OUT_MASKED
            self.outputFiles.append(outputFile)

            if self.outputGff:
                # Also merging gff3 files
                gff3File = inputDir + '/' + inputFileWithoutDir + RM_INDEX_OUT_GFF
                self.gff3files.append(gff3File)

            for index in RM_INDEX_OUTFILES_UNWANTED:
                otherOutputFile = inputDir + '/' + inputFileWithoutDir + index
                self.toDeleteFiles.append(otherOutputFile)

            countInput += 1
            stdoutFile = self.outputDir + '/' + str(countInput) + EXT_STDOUT
            stderrFile = self.outputDir + '/' + str(countInput) + EXT_STDERR
            command = ''

            # RepeatMasker -lib /ngsprod/gsap/genomes/mipsREdat_9.3p_ALL.fasta
            # -gff -x -nolow -dir /ngsprod/gsap/RunPrograms/RepeatMasker/PHIv2.1
            # /ngsprod/gsap/genomes/ZmChr1v2.fas
            command = ' '.join([
                BIN_SET_PIPEFAIL, 'cd ' + self.interimOutputDir, '&&',
                LSF_BIN_RM, self.analysisParameters, '-dir', inputDir, inputFile])

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

        # concatenation, only for masked fasta & gff3.
        # check space as the output size will be doubled temporarily.
        phi.DiskSpaceWarning.checkDiskSpace(
            self.outputDiskName, self.minFreeSpaceToWarn,
            self.minFreeSpaceToPause, self.name, EMAIL_SENDER,
            self.emails, DEFAULT_MAX_DISK_CHECK_FREQUENCY, self.logFh)

        if self.verboseLevel > 0:
            self.logger.info("Trying to concatenate the outputs into one \
file %s at %s" % (self.outputFile, phi.Utils.getTimeStampString()))

        if self.checkExistingStdoutFiles:
            # first test if the concatenated files have been
            # created before and child files are deleted.
            missingStdoutFileTotal = 0
            for f in self.stdoutFiles:
                if not os.path.exists(f):
                    missingStdoutFileTotal += 1

            if missingStdoutFileTotal == 0 and os.path.exists(self.outputFile):
                if self.verboseLevel > 0:
                    self.logger.info("Looks like concatenation was done \
before. All stdout files look ok. The output file %s is present. No rerun at \
%s" % (self.outputFile, phi.Utils.getTimeStampString()))

                self.endTime = time.localtime()
                return

        # checking & concatenating output files (.masked fasta files)
        missingOutputFileTotal = 0
        command = 'cat'
        for f in self.outputFiles:
            if not os.path.exists(f):
                self.logger.warn("Missing one of the output files %s at %s. \
Will continue merging the rest into %s." % (f, phi.Utils.getTimeStampString(), self.outputFile))
                missingOutputFileTotal += 1
            else:
                command = command + ' ' + f
        # check if no output files are created at all
        if missingOutputFileTotal == len(self.outputFiles):
            self.logger.error("No output files found from Repeat Masker. \
Check dir: %s" % self.interimOutputDir)
            raise Exception("No output files found from Repeat Masker. \
Check dir: %s" % self.interimOutputDir)

        command = ' '.join([command, '>', self.outputFile])
        if self.verboseLevel > 0:
            self.logger.info("Running concatenate command: %s" % command)
        returnCode = os.system(command)
        if returnCode != 0:
            self.logger.error("Concatenation Failed. command %s did not \
return 0: %d" % (command, returnCode))
            raise Exception("Concatenation Failed. command %s did not \
return 0: %d" % (command, returnCode))
        
        # check if sequences were sliced into chunks - if so we need to merge back
        # if (self.chunkSize > 0): - no need to check, run through mergeSeqSlices regardless
        # for consistent FASTA format, as it will ONLY merge slices
        self.logger.info("Sequence slices in output fasta file, if any, \
need to be merged back: %s" % self.outputFile)
        mergedFile = mergeSeqSlices(self.outputFile, self.logger)
        if mergedFile:
            # Overwriting output file with merged result
            command = ' '.join(['mv', mergedFile, self.outputFile])
            if self.verboseLevel > 0:
                self.logger.info("Overwriting output file with new merged file: %s" % command)
            returnCode = os.system(command)
            if returnCode != 0:
                self.logger.error("Overwriting Failed. command %s did not \
return 0: %d" % (command, returnCode))
                raise Exception("Overwriting Failed. command %s did not \
return 0: %d" % (command, returnCode))

        # merging the GFF files
        if self.outputGff:
            self.logger.info("Trying to concatenate the GFF files into one \
file %s at %s" % (self.gff3file, phi.Utils.getTimeStampString()))
            fhOut = open(self.gff3file, 'w')
            
            # in case of slices
            slice_id = re.compile(r'(.+)_slice:(\d+)-(\d+)').match

            # not simple
            for f in self.gff3files:
                # if no repeats found, no output at all - possible
                if not os.path.exists(f):
                    continue

                fhIn = open(f)
                for line in fhIn:
                    if line[0] == '#':
                        # comment lines: just copy it.
                        fhOut.write(line)
                    else:
# feature lines: re-format attr column
# Target "Motif:RLG_49411|LTR_AC198968.1_11061|LTR/Gypsy|02.01.01.10.99|4577|Zea" 15534 16097
# Re-format as
# ID=RLG_49411|LTR_AC198968.1_11061|LTR/Gypsy|02.01.01.10.99|4577|Zea:15534-16097

                        row = line.split('\t')
                        seq_id = row[0]
                        start = int(row[3])
                        end = int(row[4])
                        attrs = row[8]
                        target_attr = attrs.split(' ')
                        target_id = target_attr[1].replace('"', '') # remove quotes
                        target_id = 'ID=' + target_id.replace('Motif:', '') + ':' + \
                            target_attr[2] + '-' + target_attr[3]
                        #attrs = attrs.replace('ID=', 'ID=' + refSeqId + '_')
                        #attrs = attrs.replace('Parent=', 'Parent=' + refSeqId + '_')
                        
# look for slices and correct the seq_id, start, end values (the latter being absolute values)
#CTL1_GR2HT_1ctg_v3_slice:1500001-2000000        RepeatMasker    similarity      498999  500000  16.8 ...
#CTL1_GR2HT_1ctg_v3_slice:2000001-2500000        RepeatMasker    similarity      1       1268    20.6 ... 
                        if (self.chunkSize > 0):
                            slice_match = slice_id(seq_id)
                            if (slice_match):
                                # It is a sequence slice
                                # r'(.+)_slice:(\d+)-(\d+)'
                                seq_id = slice_match.group(1)
                                start += (int(slice_match.group(2)) - 1)
                                end   += (int(slice_match.group(2)) - 1)
                                # check if start, end of the feature are within bounds of the slice
                                if (start > int(slice_match.group(3)) or end > int(slice_match.group(3))):
                                    self.logger.error("Re-computed GFF feature bounds: {0}-{1} are not \
                                        within the slice boundaries: {2}".format(start, end, row[0]))
                                    self.logger.warn("Skipping line: %s" % line.strip())
                                    continue
                        newLine = '\t'.join([
                            seq_id, row[1], row[2], str(start), str(end),
                            row[5], row[6], row[7], target_id])
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
            if self.outputGff:
                files = [self.outputFile, self.gff3file, logFile] # + self.lsfManager.stdoutFiles
            else:
                files = [self.outputFile, logFile] # + self.lsfManager.stdoutFiles

            # somehow, outlook does not honor \n for the first a few lines.
            body = "Your command (quotes might have been lost): %s\n" % self.commandLine
            body = body + "Actual total jobs: %d\n" % self.jobTotal
            body = body + "Concurrent jobs: %d\n" % self.concurrentJobTotal
            body = body + "Cluster stats:\n"
            body = body + "CPU total: %d sec.\n" % logsObj.sumCpuTime
            body = body + "CPU max: %d sec.\n" % logsObj.maxCpuTime
            body = body + "CPU min: %d sec.\n" % logsObj.minCpuTime
            body = body + "CPU ave: %s sec.\n" % logsObj.aveCpuTimeStr
            body = body + "Memory max: %d MB\n" % logsObj.maxMemory
            body = body + "Memory min: %d MB\n" % logsObj.minMemory
            body = body + "Memory ave: %s MB\n" % logsObj.aveMemoryStr
            body = body + 'For more details, please check the files under %s\n' % self.outputDir
            body = body + 'Please also read the following result file(s) \
and log file:\n' + '\n'.join(files)

            for email in self.emails:
                phi.Utils.sendEmail(EMAIL_SENDER, email, subject, body)

        if os.path.exists(self.outputFile):
            # Delete any temporary files - if final output file is written

            # Unwanted output files from RepeatMasker - this function expects an array of files
            self.deleteFiles(self.toDeleteFiles)
            # delete child input files
            self.deleteInputFiles()
            # now delete the child output files so as to save the space.
            self.deleteOutputFiles()

            # child GFF output files
            if self.outputGff and os.path.exists(self.gff3file):
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
    """run the analysis"""

    parser = argparse.ArgumentParser(
        description="Command line tool to run Repeat Masker using HPC cluster",
        epilog=textwrap.dedent('''\
Not all combinations of parameters are supported. User is responsible \
for NOT passing conflicting parameters into the script.
                                   
Example:
python %(prog)s -i guna-testing/RepeatMasker/input/ATH_EST_sample_clean_9.fa \
-lib /ngsprod/gsap/genomes/mipsREdat_9.3p_ALL.fasta \
-dir guna-testing/RepeatMasker/output/ATH_sample -o ATH_EST_sample_clean_9 \
-gff -x -l guna-testing/RepeatMasker/output/ATH_sample/log.txt \
-P gtdi-gsap-repeatmasker -j 50 -n 50 -J RMtest -ap='-nolow' \
-lp='-R "rusage[mem=6000,scr=500]"' -s -an RM
'''),
        formatter_class=CustomFormatter)

    # Mandatory arguments
    parser.add_argument(
        '-i', dest='inputFile', required=True,
        help="**REQUIRED** input fasta file, DNA. \
It does not have to be on a cluster-accessible file system.")
    parser.add_argument(
	'-dir', dest='outputDir', required=True,
	help="**REQUIRED** Writes output to this directory. \
Will create if non-existent. Must be cluster accessible.")
    library_group = parser.add_mutually_exclusive_group(required=True)
    library_group.add_argument(
        '-species', dest='species',
        help="Specify the species or clade of the input sequence. The species name \
must be a valid NCBI Taxonomy Database species name and be contained in the \
RepeatMasker repeat database. Some commonly used species are: arabidopsis, rice, \
wheat, maize, human, mouse, mammal, carnivore, rat, cow, dog, drosophila and elegans. \
See repeatmasker -help for more info. **REQUIRED** unless -lib option is used.")
    library_group.add_argument(
        '-lib', dest='dbFile',
        help="Repeat Masker library file, DNA. Allows use of a custom library. \
**REQUIRED** unless -species option is used.")

    # Optional arguments
    parser.add_argument(
        '-o', dest='outputFilePrefix',
        help="Output files created in outputDir will start with \
this prefix name. Default: input file name prefix.")
    parser.add_argument(
        '-c', dest='chunkSize', type=int, default=0,
        help="sequence chunk size in bp. Input sequences will be \
sliced into smaller chunks of this size, if option given. Particularly \
useful for large genomes so smaller sequences can be processed in parallel \
(Ex: -c 500000 for 500 kb chunks). By default, sequences are not sliced up.")
    parser.add_argument(
        '-gff', dest='outputGff', action='store_true',
        help="Creates an additional Gene Feature Finding format output.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-x', dest='maskWithX', action='store_true',
        help="Returns repetitive regions masked with Xs rather than Ns.")
    group.add_argument(
        '-xsmall', dest='maskWithSmall', action='store_true',
        help="Returns repetitive regions in lowercase (rest capitals) rather than masked.")
    parser.add_argument(
        '-an', dest='analysisName', default=RM_ANALYSIS_NAME,
        help="analysis name. Used in email report.")
    parser.add_argument(
        '-ap', dest='analysisParameters', default=RM_ANALYSIS_PARAMETERS,
        help="RepeatMasker specific analysis parameters. Provide as one \
string within quotes with equal-to sign. Ex: -ap='-nolow'.")
    parser.add_argument(
        '-lp', dest='lsfParameters', default=LSF_RESOURCES_RM,
        help='''lsf parameters. Provide as one string within quotes with \
equal-to sign. Ex: -lp='-R "rusage[mem=6000,scr=500]"'. For special options \
only such as resources, not for queue, project name, job name etc. \
Do not use it if you are not sure.''')
    parser.add_argument(
        '-ld', dest='lsfDelayTime', type=float, default=LSF_DELAY_TIME_RM,
        help="LSF delay time in seconds. Wait this long before LSF job \
submission and after LSF job completion.")
    parser.add_argument('-q', dest='queue', default=LSF_DEFAULT_QUEUE, help="lsf queue name")
    parser.add_argument(
        '-P', dest='projectName',
        default=LSF_RM_PROJECT_NAME, help="lsf project name")
    parser.add_argument('-J', dest='jobName', default=RM_ANALYSIS_NAME, help="lsf job name")
    parser.add_argument(
        '-j', dest='maxLsfJob', type=int, default=LSF_MAX_JOB_TOTAL,
        help="max LSF jobs.")
    parser.add_argument(
        '-n', dest='maxConcurrentLsfJob', type=int,
        default=LSF_MAX_CONCURRENT_JOB_TOTAL, help="max concurrent LSF jobs.")
    parser.add_argument(
        '-l', dest='logFile',
        help="log file (not for LSF jobs), full path. Default: stderr. \
It does not have to be on a cluster-accessible file system.")
    parser.add_argument(
        '-s', '--needEmail', action='store_true',
        help="Flag set to send email upon being successfully done. \
Default: no emails. Email will always be sent out in cases of errors.")
    parser.add_argument(
        '-e', dest='emails',
        help="email addresses. Delimited by space if multiple. \
eg: -e='guna.gurazada@pioneer.com another@pioneer.com'. Default: submitter")
    parser.add_argument(
        '-z', dest='bsubInterval', type=float, default=1.0,
        help="wait time in seconds, between bsub commands. \
Can be a float, but must be less than 10.")
    parser.add_argument(
        '-r', '--requeueable', action='store_true',
        help="Flag set to requeue jobs. Default: not. If there is still a \
problem after 2 tries, the status will become PSUSP, you need to manually \
kill both the process and LSF jobs. However, there are other reasons that \
a status become PSUSP. Check LSF stdout and stderr to make sure.")
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
        '-v', dest='verboseLevel', type=int, choices=[0, 1, 2],
        default=1, help="verbose level. 0: quiet; 1: important messages; \
2: all messages. Errors will always be spit out regardless the setting.")
    parser.add_argument(
        '--version', action=CustomAction, program=LSF_BIN_RM,
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
    outputDir = args.outputDir
    species = args.species
    dbFile = args.dbFile

    if os.path.isfile(outputDir):
        outputDir = os.path.dirname(outputDir)

    outputFilePrefix = args.outputFilePrefix
    chunkSize = args.chunkSize
    outputGff = args.outputGff
    maskWithX = args.maskWithX
    maskWithSmall = args.maskWithSmall
    analysisName = args.analysisName
    analysisParameters = args.analysisParameters
    lsfParameters = args.lsfParameters
    lsfDelayTime = args.lsfDelayTime
    queue = args.queue
    projectName = args.projectName
    jobName = args.jobName
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

    # check on a few options before submitting in the try block
    if len(projectName) < 8:
        sys.stderr.write("project name: %s is too short or in wrong format.\n" % projectName)
        logFh.write("project name: %s is too short or in wrong format.\n" % projectName)
        exit(1)
    if bsubInterval > maxBsubInterval:
        sys.stderr.write("bsubInterval %f is greater than max \
%f.\n" % (bsubInterval, maxBsubInterval))
        logFh.write("bsubInterval %f is greater than max %f.\n" % (bsubInterval, maxBsubInterval))
        exit(1)
    if verboseLevel > 0:
        logFh.write("The command is as below (quotes might have been removed \
when restoring the command):\n%s\n" % ' '.join(sys.argv))

    # now comes the time-consuming part. Use try to catch the errors and report errors.
    # there are too many types of errors and impossible to check individually
    #(eg: network, space, hard drive failure). Just trace all
    try:
        # now create the object. If you are writing your own script,
        # you need to import it and use the full name:
        # phi.Analyses.RepeatMasker.RM_LSF_manyconc.RM_LSF_manyconc()
        analysisObj = RM_LSF_manyconc(
            analysisName, analysisParameters, inputFile, species, dbFile, outputFilePrefix,
            outputDir, chunkSize, outputGff, maskWithX, maskWithSmall, outputDir,
            childInputFileDir, logFh, needEmail, emails, verboseLevel, queue,
            projectName, jobName, lsfParameters, maxLsfJob,
            maxConcurrentLsfJob, checkExistingStdoutFiles, requeueable)
        analysisObj.minFreeSpaceToWarn = minFreeSpaceToWarn
        analysisObj.minFreeSpaceToPause = minFreeSpaceToPause
        analysisObj.lsfDelayTime = lsfDelayTime
        analysisObj.lsfManager.bsubInterval = bsubInterval

        analysisObj.commandLine = ' '.join(sys.argv)
        analysisObj.run()

    except:
        import traceback
        stackTrace = traceback.format_exc() + "\nSome jobs might be still \
running and you need manually bkill them\n"
        # log the error
        logFh.write(stackTrace)

        # email the users. For testing purpose, disable emailing.
        for email in emails:
            phi.Utils.sendEmail(EMAIL_SENDER, email, 'Error in ' + analysisName, stackTrace)

        raise
