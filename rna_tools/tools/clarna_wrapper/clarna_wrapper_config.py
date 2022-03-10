PATH = "/opt/ClaRNA"

try:
    from clarna_wrapper_config_local import PATH
except ImportError:
    pass