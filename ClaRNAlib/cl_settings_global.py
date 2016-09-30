import os
CLARNA_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.dirname(os.path.abspath(__file__))+"/.cache"
DATA_DIR = os.path.dirname(os.path.abspath(__file__))+"/gc-data"

PYTHON_CMD = "/usr/bin/python"
if os.path.isfile("/opt/local/bin/python"):
    PYTHON_CMD = "/opt/local/bin/python"


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# #####  You MUST set these environmental   #####
# #####  varables to be able to use these   #####
# #####  optional 3rd party programs.       #####
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


RNAVIEW_HOME="/home/mboni/bin/RNAVIEW"
# RNAVIEW_HOME="/Users/tomek/work/genesilico.programs/RNAVIEW"
# if os.path.isdir("/home/twalen/RNAVIEW"):
#     RNAVIEW_HOME="/home/twalen/RNAVIEW"
# elif os.path.isdir("/home/twalen/prgs/RNAVIEW"):
#     RNAVIEW_HOME="/home/twalen/prgs/RNAVIEW"

MCA_HOME="/home/wdawson/sw/src/mc_annotate"
# MCA_HOME="/Users/tomek/work/genesilico.programs/MC-Annotate"
# if os.path.isdir("/home/twalen/MC-Annotate"):
#     MCA_HOME="/home/twalen/MC-Annotate"
# elif os.path.isdir("/home/twalen/prgs/MC-Annotate"):
#     MCA_HOME="/home/twalen/prgs/MC-Annotate"

FR3D_HOME="/home/wdawson/sw/src/FR3D-dev"
# FR3D_HOME="/Users/tomek/work/genesilico.programs/FR3D"
# if os.path.isdir("/home/twalen/prgs/FR3D"):
#     FR3D_HOME="/home/twalen/prgs/FR3D"


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODERNA_STACKINGS_HOME = os.path.dirname(os.path.abspath(__file__))
