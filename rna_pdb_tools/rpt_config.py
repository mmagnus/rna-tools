SIMRNA_DATA_PATH = None
RCHIE_PATH = None
QRNAS_PATH = None
VARNA_PATH = None
RFAM_DB_PATH = None # path to Rfam.cm
CONTEXTFOLD_PATH = None
CPUS_CLUSTER = 1000
DIFF_TOOL="diff"
RNA_ROSETTA_RUN_ROOT_DIR_MODELING = None
QRNAS_PATH = ''
try:
    from rpt_config_local import *
except ImportError:
    pass

