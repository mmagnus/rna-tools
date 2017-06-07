BINARY_PATH = ""
try:
    from rna_x3dna_config_local import BINARY_PATH
except ImportError:
    pass

if not BINARY_PATH:
    raise Exception('Set up BINARY_PATH in rna_x3dna_config_local.py, .e.g "/Users/magnus/work/opt/x3dna/x3dna-dssr"')
