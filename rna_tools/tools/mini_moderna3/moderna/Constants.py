#!/usr/bin/env python
#
# Constants.py
#
# Contains all constant values and parameters.
# 
# http://iimcb.genesilico.pl/moderna/ 
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "rother.magdalena.gmail.com"
__status__ = "Production"

import os.path, os, math, re
import Bio
from rna_tools.tools.mini_moderna3.moderna.analyze.Constants import *
from rna_tools.tools.mini_moderna3.moderna.builder.Constants import *

#######################   path settings  ############################

MODULE_PATH, MODULE_NAME = os.path.split(__file__)
if not MODULE_PATH:
    MODULE_PATH = os.getcwd()
MODULE_PATH = re.sub('\A.*library.zip', '', MODULE_PATH) # for py2exe
if not MODULE_PATH:
    # for py2exe
    DATA_PATH = 'data' + os.sep
else:
    DATA_PATH = (MODULE_PATH + os.sep + 'data' + os.sep)


#######################   SequenceAlignment   ######################

UNKNOWN_RESIDUE_SHORT = 'X' # few-letter abbreviation
UNKNOWN_RESIDUE_ONELETTER = '.' # one-letter abbreviation
ANY_RESIDUE = UNKNOWN_RESIDUE_ONELETTER

MISSING_RESIDUE = '_' 
# Important distinction from UNKNOWN_RESIDUE_ONELETTER,
# in order to inform user that there *is* something known about that resideu
RESIDUE_WITHOUT_ONE_LETTER_ABBREV = 'x'
MODIFICATION_NAMES_TABLE_PATH = DATA_PATH + 'modification_names_table'
MODIFICATION_TOPOLOGY_FILE = DATA_PATH + 'modification_topologies.txt'
STANDARD_BASES = ['A', 'C', 'G', 'U']

###########################   Residue   ###########################

B_FACTOR_COPY = 10.00
B_FACTOR_REMOVE_MODIF = 15.00
B_FACTOR_ADD_MODIF = 20.00
B_FACTOR_EXCHANGE = 30.00
B_FACTOR_FRAGMENT = 60.00

#########################  ModernaStructure   ######################

# backbone atoms are used for superimposition  
BACKBONE_ATOMS = ['P', "O5'", "C5'", "C4'", "C3'", "O3'"]

# list of atoms to keep when removing modifications
BACKBONE_RIBOSE_ATOMS = ['P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", \
    "C3'", "O3'", "C2'", "O2'", "C1'"]
BACKBONE_RIBOSE_ATOMS_WITHOUT_O2 = ['P', 'OP1', 'OP2', "O5'", "C5'", "C4'", \
    "O4'", "C3'", "O3'", "C2'", "C1'"]

BASE_PATH = DATA_PATH + 'standard_bases/'
BASE_PAIR_PATH = '/home/krother/repos/ModeRNA/moderna/data/base_pairs/'
ADDING_MODIFICATION_RULES_PATH = DATA_PATH + 'modifications'
MODIFICATION_FRAGMENTS_PATH = DATA_PATH + 'modification_fragments/'

AA_ATOMS = ['N', 'CA', 'C', 'O']

############################  RNAModel   #########################

SINGLE_STRAND = DATA_PATH + 'single_strand.pdb'
HELIX_PATH = SINGLE_STRAND

############################ SECONDARY STRUCTURE   #########################

HELIX = DATA_PATH + 'helix.pdb'
HELIX_SUPERPOSITION = ['P', "O5'", "C5'", "C4'",  "C3'",  "O3'", 'N*']
SINGLE_PAIR = DATA_PATH +'pair_fragment.pdb'
PAIR_SUPERPOSITION = ['P', "O5'", "C5'", "C4'", "C3'", "O3'", \
    "C2'", "C1'", "O4'", 'N*']
PAIR_PURINE_SUPERPOSITION = ['C4',  'N9', 'C8']
PAIR_PYRIMIDINE_SUPERPOSITION = ['C2', 'N1', 'C6']

########################   ModernaFragment   ########################

# older
LIR_SUPERPOSITION5 = ("O3'", "C3'", "C4'", "C1'", "N*") # 04.02.2010
LIR_SUPERPOSITION3 = ("C5'", "C4'", "C3'", "C1'", "N*", "O5'") # 04.02.2010
#LIR_SUPERPOSITION5 = [ "C4'", "C3'","O3'"] # changed 04.02.2010
#LIR_SUPERPOSITION3 = ["P","O5'","C5'"] # changed 04.02.2010

############################   LIR   #############################

# stem region is a region containing continious base paring in length 
# of 4 to 11 canonical Watson-Crick base pairs.  
# stem regions with length within range
# [MINIMAL_STEM_LENGTH,MAXIMAL_STEM_LENGTH]
# won't be taken into consideration during finding fragments 
MINIMAL_STEM_LENGTH = 4 
MAXIMAL_STEM_LENGTH = 11

# [MINIMAL_LOOP_LENGTH, MAXIMAL_LOOP_LENGTH] - range of fragments lenght
MINIMAL_FRAGMENT_LENGTH = 2
MAXIMAL_FRAGMENT_LENGTH = 15

###########################   LIRdb    ##############################

MAX_FRAGMENT_LENGTH = 20
#MAX_FRAGMENT_LENGTH = 100
PATH_TO_LIR_LIST = MODULE_PATH + os.sep + 'data' + os.sep + 'rnaDB05_list.txt'

##########################   SerarchLIR   ###########################

LIR_DATABASE_PATH = MODULE_PATH + os.sep + 'data' + os.sep + 'LIR_fragments.lib'
LIR_DIRECTORY_PATH = MODULE_PATH + os.sep + 'data' + os.sep + 'rnaDB05' + os.sep
PATH_TO_LIR_STRUCTURES = LIR_DIRECTORY_PATH #TODO: cleanup duplicate name

MAX_DIST_STEM = 5.75 
NUMBER_OF_FRAGMENT_CANDIDATES = 50

# MM: check whether it is uded somewhere
# walues used for parametrization finall loopHit rank: 

#########################   CheckPdb   ###########################

PHOSPHORYLATED_NUCLEOTIDES = ['AMP', 'ADP', 'ATP', 'CMP', 'CDP', 'CTP', \
    'GMP', 'GDP', 'GTP', 'UMP',  'UDP', 'UTP']


PI2 = math.pi*2
