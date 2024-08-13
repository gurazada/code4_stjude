# by Guna

"""constants for repeatmasker"""

import os

# Alex found the problem - user perl ENV potentially conflicts with RepeatMasker tool
# user perl env can be triggerred by: source /mnt/common/profiles/ngsdev/ngsdev.env
# Line 138 and 139 in /mnt/common/profiles/ngsdev/ngsdev.env
# export PERLLIB=/mnt/common/rh6/ngsdev/lib/perllib:${PERLLIB}
# export PERL5LIB=/mnt/common/rh6/ngsdev/lib/perllib:${PERL5LIB}
# so remove them!!
if 'PERL5LIB' in os.environ:
    del os.environ['PERL5LIB']
if 'PERLLIB' in os.environ:
    del os.environ['PERLLIB']


# this is the path on the cluster, not on the submission host.
LSF_BIN_DIR = '/gsap/tools/RepeatMasker'
LSF_BIN_RM = LSF_BIN_DIR + '/RepeatMasker'


# Max Jobs to run => Input file will be split into this number of files
LSF_MAX_JOB_TOTAL = 300   # default is only 300.

# max memory usage noted is 5612MB - Wengang's PHIv2.1 tests
LSF_RESOURCES_RM = '-R "rusage[mem=7000,scr=500]"'
LSF_DELAY_TIME_RM = 120 # in seconds.
LSF_RM_PROJECT_NAME = 'gtdi-gsap-repeatmasker'

# this is for local, on submission host.
BIN_DIR = LSF_BIN_DIR
BIN_RM = BIN_DIR + '/RepeatMasker'

RM_ANALYSIS_PARAMETERS = '-nolow'


# RM output files - that should be deleted
RM_INDEX_OUTFILES_UNWANTED = ['.alert', '.cat', '.cat.gz', '.ori.out', '.out', '.tbl']

# RM output files - that need to be kept
RM_INDEX_OUT_MASKED = '.masked'
RM_INDEX_OUT_GFF = '.out.gff'


RM_INDEX_FASTA = ['.fasta', '.fas', '.fa']          # fasta extn

RM_INDEX_GFF = '.gff'


RM_GFF_FIELD2 = 'RepeatMasker'
RM_ANALYSIS_NAME = 'repeatmasker'
