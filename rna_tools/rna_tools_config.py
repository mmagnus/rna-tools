import configparser
import os

SIMRNA_DATA_PATH = None
RCHIE_PATH = None
QRNAS_PATH = None
QRNAS_CONFIG_PATH = None
VARNA_PATH = None
VARNA_JAR_NAME = None
RFAM_DB_PATH = None  # path to Rfam.cm
CONTEXTFOLD_PATH = None
IPKNOT_PATH = "/usr/local/opt/bin/ipknot"
CPUS_CLUSTER = 1000
DIFF_TOOL = "diff"

RNA_ROSETTA_RUN_ROOT_DIR_MODELING = ""
RNA_ROSETTA_NSTRUC = 10000

EASY_CAT_PATH = ""
RNASTRUCTURE_PATH = ""
ENTRNA_PATH = ''
#QRNAS_PATH = os.getenv('QRNAS_PATH', PATH + '/opt/qrnas/')
LOG_DIRECTORY = ''

RASP_PATH = ""       # "/Users/magnus/work/papers/mqaprna/mqaprna_env/mqapRNA/opt/rasp-fd-1.0/"
RASP_PATH_DIR = ""   # /Users/magnus/work/papers/mqaprna/mqaprna_env/mqapRNA/opt/rasp-fd-1.0"

LOG_DIRECTORY = ""   # "/tmp/"
TMP_PATH = ""
BIN_PATH = ""

WRAPPERS_PATH = ""   # "/Users/magnus/work/src/rna-tools/rna_tools/tools/mq/"

dfire_PATH = ""      # "/Users/magnus/work/papers/mqaprna/mqaprna_env/mqapRNA/opt/dfire_rna/"

RNA3DCNN_PATH=""     # "/Users/magnus/work/papers/mqaprna/mqaprna_env/mqapRNA/opt/RNA3DCNN"
PYTHON3_PATH=""      # "/Users/magnus/miniconda2/envs/py37/bin/python"
FARNA_PATH=""        # "/Users/magnus/work/opt/rosetta_bin_mac_2018.33.60351_bundle/main/source/build/src/release/macos/10.13/64/x86/clang/9.0/static/rna_minimize.static.macosclangrelease"
FARNA_DB_PATH=""     # "/Users/magnus/work/opt/rosetta_bin_mac_2018.33.60351_bundle/main/database"
FARNA_LORES=""       # "/scoring/weights/rna/denovo/rna_lores.wts"

FARFAR2_PATH=""
FARFAR2_DB_PATH=""
FARFAR2_LORES=""

QRNA_PATH=""         # "/Users/magnus/work/opt/qrnas/"
QRNA_CONFIG_PATH=""  # "/Users/magnus/work/opt/qrnas/"
PHENIX_BIN_PATH=""   # /Applications/phenix-1.18.2-3874/build/bin"
baRNAba_PATH=""      # /Users/magnus/work/opt/barnaba/barnaba_201128"
RNAscore_PATH=""
WRAPPERS_PATH=""     # /Users/magnus/work/src/rna-tools/rna_tools/tools/mq
GRMLIB_PATH=""       # /usr/local/gromacs/share/gromacs/top/
DB_AA_PATH=""        # /Users/magnus/work/papers/mqaprna/mqaprna_env/db/RNA_aa_full
DB_5PT_PATH=""       # /Users/magnus/work/papers/mqaprna/mqaprna_env/db/RNA_5pt_full_sc1'
GROMACS_LD_PATH=""
GROMACS_PATH_BIN=""
SIMRNA_PATH=""
SIMRNA_DATA_PATH=""

X3DNA=""
X3DNA_FP=""

baRNAba_data_PATH = ''

# https://docs.python.org/3/library/configparser.html
# based on config
## def is_int(s):
##     try:
##         int(s)
##         return True
##     except ValueError:
##         return False
## config = configparser.ConfigParser()
## 
## conf = config.read(user_path + '/.rna_tools.ini')
## if conf:
##     for var in config['rna-tools'].keys():
##         value = config['rna-tools'][var]
##         if is_int(value):
##             globals()[var.upper()] = int(value)
##         else:
##             globals()[var.upper()] = value
user_path = os.path.expanduser("~")
try:
    exec(open(user_path + '/.rna_tools.py').read())  # python3
except: # FileNotFoundError: noooot perfect! 
    pass
# for cruel testing
# print(DIFF_TOOL, '\n', VARNA_PATH) # globals()
