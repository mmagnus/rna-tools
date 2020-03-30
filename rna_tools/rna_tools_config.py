import configparser
import os

SIMRNA_DATA_PATH = None
RCHIE_PATH = None
QRNAS_PATH = None
VARNA_PATH = None
VARNA_JAR_NAME = None
RFAM_DB_PATH = None  # path to Rfam.cm
CONTEXTFOLD_PATH = None
CPUS_CLUSTER = 1000
DIFF_TOOL = "diff"

RNA_ROSETTA_RUN_ROOT_DIR_MODELING = ""
RNA_ROSETTA_NSTRUC = 10000

EASY_CAT_PATH = ""
RNASTRUCTURE_PATH = ""
ENTRNA_PATH = ''
#QRNAS_PATH = os.getenv('QRNAS_PATH', PATH + '/opt/qrnas/')

user_path = os.path.expanduser("~")

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
try:
    exec(open(user_path + '/.rna_tools.py').read())  # python3
except: # FileNotFoundError: noooot perfect! 
    pass
# for cruel testing
# print(DIFF_TOOL, '\n', VARNA_PATH) # globals()
