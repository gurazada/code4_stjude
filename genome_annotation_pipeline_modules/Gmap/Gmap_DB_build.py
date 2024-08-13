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
gmap_db_build.py -d /anno/gsap/benchmark/TAIR10_chr_all.fa -D \
/ngsprod/gsap/GMAPDB -db TAIR10_chr_all -l log.txt -ap='-k 14 -s chrom' \
-s -an GMAP_BUILD
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
from phi.Analyses.Constants import DEFAULT_MIN_FREE_SPACE_TO_WARN, \
    DEFAULT_MIN_FREE_SPACE_TO_PAUSE, EMAIL_SENDER, \
    DEFAULT_MAX_DISK_CHECK_FREQUENCY
import phi.Utils
import phi.Analyses.Local.Analysis2
import phi.DiskSpaceWarning2
import phi.Logger

from .Gmap_Constants import GMAP_BUILD_ANALYSIS_NAME, \
    GMAP_BUILD_ANALYSIS_PARAMETERS, LSF_DIR_GMAPDB, \
    GMAP_INDEX_FASTA, GMAP_INDEX_FAS, GMAP_INDEX_FA, GMAP_INDEX_DBFILES, \
    LSF_BIN_GMAP_BUILD_V1, LSF_BIN_GMAP_BUILD_V2
from phi.Parse import cleanFastaNSplit

# Setting Global variables for GMAP program paths - default V1
LSF_BIN_GMAP_BUILD = LSF_BIN_GMAP_BUILD_V1

class Local(phi.Analyses.Local.Analysis2.Analysis):
    """GMAP DB build module that runs gmap_build locally (not via LSF)."""

    def __init__(
            self,
            name=GMAP_BUILD_ANALYSIS_NAME,
            analysisParameters=GMAP_BUILD_ANALYSIS_PARAMETERS,
            dbName='',
            dbDestDir=LSF_DIR_GMAPDB,
            dbFile='',
            needEmail=True,
            emails=[phi.Utils.getInteralEmailAddress()],
            logger=None
        ):

        phi.Analyses.Local.Analysis2.Analysis.__init__(
            self, name, LSF_BIN_GMAP_BUILD, analysisParameters, '',
            dbFile, [[]], [[]], '', dbDestDir, '', needEmail, emails,
            None, logger)

        # gmap_build specific properties
        self.dbName = dbName
        self.dbDestDir = dbDestDir
        self.dbFile = dbFile
        self.toDeleteFiles = [] # Add any temporary files to this array, that are SURE to go.

    def run(self):
        # check on db file - required
        if not self.dbFile or not os.path.exists(self.dbFile):
            self.logger.error("DB fasta file not set or does not exist: %s" % self.dbFile)
            raise Exception, "DB fasta file not set or does not exist: %s" % self.dbFile
        self.dbFile = os.path.abspath(self.dbFile)

        # -D option for GMAP
        dbDestDir = self.dbDestDir
        if not os.path.exists(dbDestDir):
            self.logger.error("GMAP DB dest dir does not exist: %s" % dbDestDir)
            raise Exception, "GMAP DB dest dir does not exist: %s" % dbDestDir
        self.dbDestDir = os.path.abspath(dbDestDir)

        if not os.access(self.dbDestDir, os.W_OK) or not os.access(self.dbDestDir, os.X_OK):
            self.logger.error("GMAP DB dest dir %s is either not writable or \
executable. I need both to be able to write to it" % self.dbDestDir)
            raise Exception, "GMAP DB dest dir %s is either not writable or \
executable. I need both to be able to write to it" % self.dbDestDir


        # -d option for GMAP (-db here)
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
            # check space:
            phi.DiskSpaceWarning2.checkDiskSpace(
                self.dbDestDir, self.minFreeSpaceToWarn,
                self.minFreeSpaceToPause, self.name, EMAIL_SENDER, self.emails,
                DEFAULT_MAX_DISK_CHECK_FREQUENCY, self.logger)

            self.logger.info("GMAP DB does not exist or missing index files. \
Need to create first at: %s" % dbIndexFilesDir)
            self.logger.info("Parsing DB fasta file: %s..." % self.dbFile)
            # Parse DB fasta file to clean; No splitting -
            # returns an array still. Empty if invalid sequences
	    seqNameFile = self.dbDestDir + '/' + self.dbName + '_origSeqNames'
            newDbFiles = cleanFastaNSplit(
                self.dbFile, 1, seqNameFile, 0, self.dbDestDir, self.dbName, self.logger)
            if len(newDbFiles) == 0:
                self.logger.error("DB fasta file has no valid sequences: %s" % self.dbFile)
                raise Exception, "DB fasta file has no valid sequences: %s" % self.dbFile
            # Run GMAP build
            # gmap_build -d TAIR10_chr_all -D /ngsprod/gsap/GMAPDB
            #   /anno/gsap/benchmark/TAIR10_chr_all.fa
            command = ' '.join([
                self.bin, '-d', self.dbName, '-D', self.dbDestDir,
                self.analysisParameters, newDbFiles[0]
            ])
            self.logger.info("Running GMAP DB build command: %s" % command)
            returnCode = os.system(command)
            if returnCode != 0:
                self.logger.error(
                    "DB build Failed. command %s did not return 0: %d" % (command, returnCode)
                )
                raise Exception, "DB build Failed. command %s did not \
return 0: %d" % (command, returnCode)
            self.logger.info("GMAP DB build done: %s" % dbIndexFilesDir)
            self.toDeleteFiles.append(newDbFiles[0])


        self.endTime = time.ctime()
        self.logger.info("Analysis All well done")

        self.postProcess()

    def postProcess(self):
        """post process"""

        if self.needEmail:
            subject = 'Analysis %s is well done at %s' % (self.name, self.endTime)
            body = "Your command (quotes might have been lost): %s\n" % self.commandLine
            body = body + "Your output files at: %s\n" % (self.dbDestDir + '/' + self.dbName)
            for email in self.emails:
                phi.Utils.sendEmail(EMAIL_SENDER, email, subject, body)

        # Delete any temporary files (Ex: Newly created DB file) -
        # this function expects an array of files
        self.deleteFiles(self.toDeleteFiles)
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
                print(cmdline(prog + ' --version | head -4'))
        else:
	    print("Analysis program invoked: %s" %(program))
	    print(cmdline(program + ' --version | head -4'))
        #message = cmdline(program + ' --version')
        #formatter = parser._get_formatter()
        #formatter.add_text(message)
        #parser._print_message(formatter.format_help(), argparse._sys.stdout)
        parser.exit()

def run():
    """run the analysis"""

    parser = argparse.ArgumentParser(
        description="Command line tool to build GMAP DB files on local machine.",
        epilog=textwrap.dedent('''\
        Not all combinations of parameters are supported. \
        User is responsible for NOT passing conflicting parameters into the script.
                                   
        Example:
        python %(prog)s -d /anno/gsap/benchmark/TAIR10_chr_all.fa -D /ngsprod/gsap/GMAPDB -db TAIR10_chr_all -l log.txt -ap='-k 14 -s chrom' -s -an GMAP_BUILD
        '''),
        formatter_class=CustomFormatter)

    # Mandatory arguments
    parser.add_argument(
        '-d', dest='dbFile', required=True, help="**REQUIRED** database file, \
DNA. If GMAP database doesn't exist, the program will try to build (-D, -db \
options). The filename part should be xxx.fa. xxx part should be alphanumeric.")

    # Optional arguments
    parser.add_argument('-D', dest='dbDestDir', default=LSF_DIR_GMAPDB, help="destination \
dir for GMAP database. If provided, dir should exist. GMAP DB for the dbFile \
will be built here if missing index files.")
    parser.add_argument('-db', dest='dbName', help="GMAP database name \
(sub-directory within dbDestDir). Will be prefix of dbFile if not given.")
    parser.add_argument(
        '-an', dest='analysisName',
        default=GMAP_BUILD_ANALYSIS_NAME, help="analysis name. Used in email report.")
    parser.add_argument(
        '-ap', dest='analysisParameters', default=GMAP_BUILD_ANALYSIS_PARAMETERS,
        help="GMAP BUILD specific analysis parameters. Provide as one string within \
quotes with equal-to sign. Ex: -ap='-k 14 -s chrom'. \
Does not conflict with other program-specific options.")
    parser.add_argument('-l', dest='logFile', help="log file, full path. Default: stderr.")
    parser.add_argument(
        '-s', '--needEmail', action='store_true',
        help="Flag set to send email upon being successfully done. Default: \
no emails. Email will always be sent out in cases of errors.")
    parser.add_argument(
        '-e', dest='emails', help="email addresses. Delimited by space if multiple. \
eg: -e='guna.gurazada@pioneer.com another@pioneer.com'. Default: submitter")
    parser.add_argument(
        '-m', dest='minFreeSpaceToWarn', type=int,
        default=DEFAULT_MIN_FREE_SPACE_TO_WARN/1000000,
        help="min free disk space to warn. Integer in MB.")
    parser.add_argument(
        '-M', dest='minFreeSpaceToPause', type=int,
        default=DEFAULT_MIN_FREE_SPACE_TO_PAUSE/1000000,
        help="min free disk space to pause. Integer in MB.")
    parser.add_argument(
        '-V', '--useVersion', dest='gmapVersion', type=int, choices=[1, 2], default=1,
        help="GMAP BUILD version to run. Use --version flag to see which GMAP versions are installed \
and set this option accordingly. Make sure to use appropriate -D dbDestDir option for the version selected.")
    parser.add_argument(
        '--version', action=CustomAction, program=[LSF_BIN_GMAP_BUILD_V1, LSF_BIN_GMAP_BUILD_V2],
        help="show analysis program's version and exit")

    args = parser.parse_args()

    # open log file first
    logFh = None
    if not args.logFile:
        logFh = sys.stderr
    else:
        logFh = open(args.logFile, 'a', 0) # no buffering.
    logger = phi.Logger.Logger(args.analysisName, logFh)

    # parsing options
    dbFile = args.dbFile
    dbDestDir = args.dbDestDir
    dbName = args.dbName
    analysisName = args.analysisName
    analysisParameters = args.analysisParameters

    needEmail = args.needEmail
    if not args.emails:
        emails = [phi.Utils.getInteralEmailAddress()]
    else:
        emails = args.emails.split(' ')
    minFreeSpaceToWarn = args.minFreeSpaceToWarn * 1000000  # input is in mega
    minFreeSpaceToPause = args.minFreeSpaceToPause * 1000000 # input is in mega
    gmapVersion = args.gmapVersion
    
    # change GMAP path based on the version provided
    # default is set for version 1 - so only change if version 2 is given 
    if gmapVersion == 2:
	global LSF_BIN_GMAP_BUILD
	LSF_BIN_GMAP_BUILD = LSF_BIN_GMAP_BUILD_V2

    # running command
    logger.info("The command is as below (quotes might have been removed \
when restoring the command):\n%s" % ' '.join(sys.argv))

    # now comes the time-consuming part.
    # Use try to catch the errors and report errors.
    # there are too many types of errors and impossible to check individually
    # (eg: network, space, hard drive failure). Just trace all
    try:
        # now create the object. If you are writing your own script, you need
        # import it and use the full name: phi.Analyses.Gmap.Gmap_DB_build.Local()
        analysisObj = Local(
            analysisName, analysisParameters, dbName, dbDestDir, dbFile, needEmail, emails, logger)
        analysisObj.minFreeSpaceToWarn = minFreeSpaceToWarn
        analysisObj.minFreeSpaceToPause = minFreeSpaceToPause
        analysisObj.commandLine = ' '.join(sys.argv)
        analysisObj.run()

    except:
        import traceback
        stackTrace = traceback.format_exc()  # + any additional info.
        # log the error
        logger.error(stackTrace)

        # email the users. For testing purpose, disable emailing.
        for email in emails:
            phi.Utils.sendEmail(EMAIL_SENDER, email, 'Error in ' + analysisName, stackTrace)

        raise



