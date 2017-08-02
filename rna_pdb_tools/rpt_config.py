SIMRNA_DATA_PATH = None
RCHIE_PATH = None
QRNAS_PATH = None
VARNA_PATH = None
RFAM_DB_PATH = None  # path to Rfam.cm
CONTEXTFOLD_PATH = None
CPUS_CLUSTER = 1000
DIFF_TOOL = "diff"

RNA_ROSETTA_RUN_ROOT_DIR_MODELING = "/home/magnus/rosetta-runs"
RNA_ROSETTA_NSTRUC = 10000

EASY_CAT_PATH = ""

import os
try:
    PATH = os.environ['RNA_PDB_TOOLS']
except KeyError:
    print ('Set up RNA_PDB_TOOLS, see Installation note')
    pass
else:
    QRNAS_PATH = os.getenv('QRNAS_PATH', PATH + '/opt/qrnas/')

try:
    from rpt_config_local import *
except ImportError:
    pass
