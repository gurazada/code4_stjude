# by Guna
"""Constant for GMAP"""

# this is the path on the cluster, not on the submission host.
# V1
LSF_BIN_DIR_V1 = '/mnt/common/rh6/cluster/tdi/bin'
LSF_BIN_GMAP_V1 = LSF_BIN_DIR_V1 + '/gmap'
LSF_BIN_GMAP_BUILD_V1 = LSF_BIN_DIR_V1 + '/gmap_build'

# V2
# need to install both LOCAL & LSF copies for V2, keeping consistent
# gmap uses LSF only & gmap_build uses LOCAL only
# gmap LOCAL copy is only used for --version flag, since that runs locally
#LOCAL_BIN_DIR_V2 = '/gsap/tools/sandbox/bin' # test install
LOCAL_BIN_DIR_V2 = '/gsap/tools/bin'
#LSF_BIN_DIR_V2 = '/mnt/common/tdi/bin' # cluster test install
LSF_BIN_DIR_V2 = '/mnt/common/tdi/gsap/bin'
LSF_BIN_GMAP_V2 = LSF_BIN_DIR_V2 + '/gmap'
LOCAL_BIN_GMAP_V2 = LOCAL_BIN_DIR_V2 + '/gmap'
LSF_BIN_GMAP_BUILD_V2 = LOCAL_BIN_DIR_V2 + '/gmap_build'

LSF_DIR_GMAPDB = '/ngsprod/gsap/GMAPDB'

# Max Jobs to run => Input file will be split into this number of files
LSF_MAX_JOB_TOTAL = 300   # default is only 300.

# max memory usage noted is 16000MB - new GMAP versions
LSF_RESOURCES_GMAP = '-R "rusage[mem=20000,scr=100]"'
LSF_DELAY_TIME_GMAP = 120 # in seconds.
LSF_GMAP_PROJECT_NAME = 'gtdi-gsap-gmap'

GMAP_ANALYSIS_PARAMETERS = '-F -n 20 -K 50000 --min-identity=0.95 --min-trimmed-coverage=0.90'
GMAP_BUILD_ANALYSIS_PARAMETERS = '-s chrom'

# one of the 8 GMAP DB index files - that SHOULD be there
GMAP_INDEX_DBFILES = [
    '.chromosome.iit',
    '.contig.iit',
    '.genomebits128',
    '.genomecomp',
    '.ref153offsets64meta',
    '.ref153offsets64strm',
    '.ref153positions',
    '.version'
]

GMAP_INDEX_FASTA = '.fasta'          # fasta extn
GMAP_INDEX_FA = '.fa'             # fasta extn
GMAP_INDEX_FAS = '.fas'            # fasta extn

GMAP_INDEX_GFF3 = '.gff3'

# Unique IDs: gff3_gene, gff3_match_est; Non-unique IDs: gff3_match_cdna
# gff3_gene is default
GMAP_OUTPUT_FORMAT_GFF3 = ['gff3_gene', 'gff3_match_cdna', 'gff3_match_est']


GMAP_GFF3_FIELD2 = 'GMAP'
GMAP_ANALYSIS_NAME = 'gmap'
GMAP_BUILD_ANALYSIS_NAME = 'gmap_build'
GMAP_TRANSCRIPT_TYPE = ['EST', 'cDNA', 'Assembly']  # Assembly is default - last value in array
