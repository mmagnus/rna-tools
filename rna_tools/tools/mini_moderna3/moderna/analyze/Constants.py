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
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

import math

#######################   RNA commons   ############################

# atoms essential for nucleotides.
# O2' not listed because of DNA
RIBOSE = ["C1'", "C2'", "C3'", "C4'", "C5'", "O3'", "O4'", "O5'"]
NUCLEOTIDE_ATOMS = {
    "A":set(RIBOSE + ["N1", "C2", "N3", "C4", "C5", "C6", \
            "N6", "N7", "C8", "N9",]),
    "G":set(RIBOSE + ["N1", "C2", "N2", "N3", "C4", "C5", \
            "C6", "O6", "N7", "C8", "N9"]),
    "C":set(RIBOSE + ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]),
    "U":set(RIBOSE + ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]),
    "T":set(RIBOSE + ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C5M", "C6"])
    }

PHOSPHATE_GROUP = set(['P', 'OP1', 'OP2', 'OP3'])

STANDARD_BASES = list(NUCLEOTIDE_ATOMS.keys())

AMINO = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
         'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',            
         'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
         'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y',
         'MSE':'m',
         }

# List is incomplete.
HETERO_GROUPS = ['ZN', 'NA', 'CO', 'FMT', 'SO4', \
    'PO4', 'SUC', 'SR', 'MN', 'CL']


#######################   HBondCalculator   ##########################

# MAXIMAL DISTANCE BETWEEN C1 OF RESIDUES
MAX_C1_DIST = 15.0

# length of covalent bond to hydrogen
H_COVALENT_BOND = 1.1
H_ANGLE_ONE = 120
H_ANGLE_TWO = 120

# parameters for hydrogen sampling
H_GENERATE_STEP = 30
H_GENERATE_ANGLE = 60.0
H_GENERATE_RADIUS = math.sin(math.pi*H_GENERATE_ANGLE/180) * H_COVALENT_BOND
H_GENERATE_COS = math.cos(math.pi*H_GENERATE_ANGLE/180) * H_COVALENT_BOND

# hbond boundaries
MAXD1 = 2.8 # Suehnel paper gives 2.5 but we tried the more tolerant 2.8
MAXD2 = 3.9 # Suehnel paper gives 3.9

# donor/acceptor atoms taken from RNAview
DONORS = {
    'G':["N1","C8","N2","O2'",],
    'A':["N6","C8","C2","O2'",],
    'I':["N1","C8","C2","O2'",],
    'C':["N4","C5","C6","O2'"],
    'U':["N3","C5","C6","O2'"],
    'P':["N1","N3","C6","O2'"],
    'T':["N3","C5","C6"],
    }

ACCEPTORS = {
    'G':["O2'","O3'","O4'","O5'","OP1","OP2","N9","N7","O6","N3","N1"],
    'A':["O2'","O3'","O4'","O5'","OP1","OP2","N9","N7","N3","N1"],
    'I':["O2'","O3'","O4'","O5'","OP1","OP2","N9","N7","O6","N3","N1"],
    'C':["O2'","O3'","O4'","O5'","OP1","OP2","N1","O2", "N3"],
    'U':["O2'","O3'","O4'","O5'","OP1","OP2","O4","O2","N1"],
    'P':["O2'","O3'","O4'","O5'","OP1","OP2","C5","O4","O2"],
    'T':["O2'","O3'","O4'","O5'","OP1","OP2","O4","O2","N1"],
    }
    
# atom names to find neighbors
PURINE_NEIGHBOR_TABLE = {
    " OP1":["P"], 
    " OP2":["P"], 
    " C2'":["O2'", "C1'", "C3'"], 
    " O2'":["C2'"], 
    " O3'":["C3'"], 
    " O4'":["C4'", "C1'"], 
    " O5'":["P", "C5'"], 
    
    " N1 ":["C6", "C2"], 
    " N2 ":["C2"], 
    " C2 ":["O2", "N1", "N2", "N3"], 
    " N3 ":["C2","C4"], 
    " C4 ":["C5", "N9", "N3"], 
    " C5 ":["C6", "C4", "N7"], 
    " C6 ":["N1", "C5", "N6", "O6"], 
    " O6 ":["C6"], 
    " N6 ":["C6"],  
    " N7 ":["C5","C8"], 
    " C8 ":["N7", "N9"], 
    " N9 ":["C8", "C4", "C1'"], 
    }
PYRIMIDINE_NEIGHBOR_TABLE = {
    " OP1":["P"], 
    " OP2":["P"], 
    " C2'":["O2'", "C1'", "C3'"], 
    " O2'":["C2'"], 
    " O3'":["C3'"], 
    " O4'":["C4'", "C1'"], 
    " O5'":["P", "C5'"], 
    
    " N1 ":["C1'","C6", "C2"], 
    " O2 ":["C2"], 
    " N2 ":["C2"], 
    " C2 ":["O2", "N1", "N2", "N3"], 
    " N3 ":["C2","C4"], 
    " O4 ":["C4"], 
    " N4 ":["C4"], 
    " C4 ":["C5", "N4", "O4", "N3"], 
    " C5 ":["C6", "C4"], 
    " C6 ":["N1", "C5"], 
    }

#######################  Base Pairs   ##########################

# for determining standard base pairs
WC_BASE_PAIRS = {
        'A':'U', 
        'G':'C', 
        'C':'G', 
        'U':'A'
        }

#######################   StackingCalculator   ##########################

# between purines and pyrimidines, the normals are reversed, because 
# the rings are reversed with respect to the helix axis. 
NORMAL_SUPPORT = {
        'C':['N1','C2','N3','C4','C5','C6'],
        'U':['N1','C2','N3','C4','C5','C6'],
        'T':['N1','C2','N3','C4','C5','C6'],
        'G':['N1','C2','C4','N3','C5','C6'],
        'A':['N1','C2','C4','N3','C5','C6'],   
        }

# reverse stacking
REVERSE = {
    '>>':'<<',
    '<<':'>>',
    '<>':'<>',
    '><':'><',
    }

ARCPI = 180.0/math.pi

#######################   ClashRecognizer   ##########################

ATOM_RADII = {'C' : 0.72, 'N' : 0.69, 'P' : 0.97, 'O' : 0.52, 'S' : 1.04}
# All radii above were taken from 
# http://geometry.molmovdb.org/files/
# libproteingeometry/data/NucProt.atom-defs.dat
# C radius corresponds to C3H0s radius
# N radius corresponds to N3H0u radius
# P radius corresponds to P4H0u radius
# O radius corresponds to O1H0u radius
# S radius corresponds to S2H0u and S2H1u radii

SEARCH_RADIUS = 2.1 # max sum of two radii        

