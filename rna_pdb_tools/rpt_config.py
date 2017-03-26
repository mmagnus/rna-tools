SIMRNA_DATA_PATH = None
RCHIE_PATH = None
QRNAS_PATH = None
RFAM_DB_PATH = None # path to Rfam.cm

try:
    from rpt_config_local import *
except ImportError:
    pass
