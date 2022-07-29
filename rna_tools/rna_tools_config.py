"""
Example::

   from rna_tools.rna_tools_config import PYMOL_PATH

"""
import configparser
import os
import pathlib
RT = os.path.dirname(str(pathlib.Path(__file__).parent.resolve())) + os.sep

SIMRNA_DATA_PATH = None
RCHIE_PATH = None
QRNAS_PATH = RT + '/opt/qrnas/'
QRNAS_CONFIG_PATH = RT + '/opt/qrnas/'
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

# MQ #
METHOD_LIST =  ['ClashScore', 'AnalyzeGeometry', 'SimRNA_0',   'RNAscore', 'eSCORE','RNAkb',
                 'RASP', 'RNAkb_all', 'RNA3DCNN', 'Dfire', 'FARNA', 'FARNA_hires', 'FARFAR2']
# default options for wrappers
WRAPPER_OPTIONS = dict([(m, []) for m in METHOD_LIST])
try:
    WRAPPER_OPTIONS['RASP'].append('all')    # all atom representation
except KeyError:
    pass
try:
    WRAPPER_OPTIONS['RNAkb'].append('5pt')  
except KeyError:
    pass
ML_MODEL_PATH = ''

RASP_PATH = RT + 'opt' + os.sep + "rasp-fd-1.0/"
RASP_PATH_DIR =  RT + 'opt' + os.sep + "rasp-fd-1.0/" 

tmp = RT + 'tmp' + os.sep
try:
    os.mkdir(tmp)
except:
    pass

LOG_DIRECTORY = tmp
TMP_PATH = tmp
SANDBOX_PATH = tmp
BIN_PATH = "/opt/homebrew/bin/"

WRAPPERS_PATH = RT + "/rna_tools/tools/mq/"

dfire_PATH = RT + "/opt/dfire_rna/"      # "/Users/magnus/work/papers/mqaprna/mqaprna_env/mqapRNA/opt/dfire_rna/"

RNA3DCNN_PATH=""     # "/Users/magnus/work/papers/mqaprna/mqaprna_env/mqapRNA/opt/RNA3DCNN"
PYTHON3_PATH=""      # "/Users/magnus/miniconda2/envs/py37/bin/python"
FARNA_PATH=""        # "/Users/magnus/work/opt/rosetta_bin_mac_2018.33.60351_bundle/main/source/build/src/release/macos/10.13/64/x86/clang/9.0/static/rna_minimize.static.macosclangrelease"
FARNA_DB_PATH=""     # "/Users/magnus/work/opt/rosetta_bin_mac_2018.33.60351_bundle/main/database"
FARNA_LORES=""       # "/scoring/weights/rna/denovo/rna_lores.wts"

opt = RT + 'opt' + os.sep
import platform
if platform.system() == "Darwin":
    FARFAR2_PATH = opt + 'ff2/rna_minimize.static.macosclangrelease'
else:
    FARFAR2_PATH = opt + 'ff2/rna_minimize.static.linuxgccrelease'
    
FARFAR2_DB_PATH = opt + 'ff2/database'
FARFAR2_LORES = ""


QRNA_PATH=""         # "/Users/magnus/work/opt/qrnas/"
QRNA_CONFIG_PATH=""  # "/Users/magnus/work/opt/qrnas/"
PHENIX_BIN_PATH=""   # /Applications/phenix-1.18.2-3874/build/bin"
baRNAba_PATH=""      # /Users/magnus/work/opt/barnaba/barnaba_201128"
RNAscore_PATH=""
GRMLIB_PATH=""       # /usr/local/gromacs/share/gromacs/top/
DB_AA_PATH=""        # /Users/magnus/work/papers/mqaprna/mqaprna_env/db/RNA_aa_full
DB_5PT_PATH=""       # /Users/magnus/work/papers/mqaprna/mqaprna_env/db/RNA_5pt_full_sc1'
GROMACS_LD_PATH=""
GROMACS_PATH_BIN=""
SIMRNA_PATH=""
SIMRNA_DATA_PATH=""
rsRNASP_PATH = RT + '/opt/rsRNASP/'


X3DNA= RT + 'opt' + os.sep + 'dssr/x3dna-dssr-64bit' #x3dna-dssr-64bit
X3DNA_FP=""

baRNAba_data_PATH = ''

PYMOL_PATH = '/usr/lib/python3/dist-packages/'
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

