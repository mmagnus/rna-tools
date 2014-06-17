PLATFORM = '64bit' # default
try:
    from py3dna_config_local import PLATFORM
except ImportError:
    pass