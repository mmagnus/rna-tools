#!/usr/bin/env python
#
# mol_parameters.py
#
# Parameters for annotating molecular structures.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, Kristian Rother"
__credits__ = ["Sabrina Hofmann"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

DEFAULT_BOND_LENGTHS = ( (1.4,2.0), (1.2, 1.4), (1.1, 1.2))

MAX_BOND_DISTANCE = 2.8
# distance ranges for detecting single bonds in 3D coordinates
SINGLE_BOND_DISTANCES = {
    'C' :{'C':(1.41,1.9),'O':(1.35,1.71),'N':(1.38,1.65),'P':(1.4,1.9),'S':(1.65,2.0),'Se':(1.5,2.5),},
    'O' :{'C':(1.35,1.71),'O':(1.4,1.6),'N':(1.4,1.6),'P':(1.56,1.74),'S':(1.7,2.0),'Se':(0.0,0.0),},
    'N' :{'C':(1.38,1.65),'O':(1.4,1.6),'N':(1.4,1.6),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    'P' :{'C':(1.4,1.9),'O':(1.56,1.74),'N':(0.0,0.0),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    'S' :{'C':(1.65,2.0),'O':(1.7,2.0),'N':(0.0,0.0),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    'Se':{'C':(1.5,2.5),'O':(0.0,0.0),'N':(0.0,0.0),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    }

# distance ranges for detecting double bonds in 3D coordinates
DOUBLE_BOND_DISTANCES = {
    #CC: 1.25
    'C' :{'C':(1.19,1.41),'O':(1.05,1.35),'N':(1.16,1.38),'P':(0.0,0.0),'S':(1.2,1.65),'Se':(1.2,1.8),},
    'O' :{'C':(1.05,1.35),'O':(1.15,1.40),'N':(1.15,1.40),'P':(1.2,1.56),'S':(1.5,1.7),'Se':(0.0,0.0),},
    'N' :{'C':(1.16,1.38),'O':(1.15,1.40),'N':(1.15,1.40),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    'P' :{'C':(0.00,0.00),'O':(1.2,1.56),'N':(0.00,0.00),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    'S' :{'C':(1.20,1.65),'O':(1.5,1.7),'N':(0.00,0.00),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    'Se':{'C':(1.40,1.80),'O':(0.0,0.00),'N':(0.00,0.00),'P':(0.0,0.0),'S':(0.0,0.0),'Se':(0.0,0.0),},
    }

TRIPLE_BOND_DISTANCES = {
    'C':{'N':(1.10,1.16)},'N':{'C':(1.10,1.16)}
    }

# electronegativities of elements
# taken from the Pauling scale, see wikipedia.
# the Rxx elements are carbon residues as they occur in KEGG
ELECTRONEGATIVITIES = {
    'C':2.55, 'O':3.44, 'N':3.04, 'H':2.20,
    'R':2.55, 'S':2.58, 'P':2.19, 'Na':0.93,
    'F':3.98, 'Cl':3.16, 'Br':2.96, 'I':2.66,
    'Cu':1.90, 'Fe':1.83, 'Mg':1.31, 'As':2.18,
    'Ca':1.00, 'K':0.82, 'Li':0.98, 'Se':2.55,
    'Mn':1.55, 'Zn':1.65, 'Mo':2.16, 'Co':1.88, 
    'Ni':1.91, 'Cr':1.66, 'V':1.63, 'Hg':2.00, 
    'W':2.36, 'Cd':1.69, 'Te':2.1, 'Pb':2.33,
    'Ag':1.93, 'Au':2.54, 'Pt':2.28, 'Si':1.90,
    'Sb':2.05, 'Bi':2.02, 'Ga':1.81, 'Sn':1.96,
    'Tc':1.9, 'B':2.04, 'Gd':1.2, 'Al':1.61,
    'Ba':0.89,'Cs':0.7,'Tl':1.62,'Pd':2.20,
    'Be':1.0,'In':1.0,'Ru':1.0,'Rh':1.0,'Eu':1.0,'La':1.0,
    'Ce':1.0,'Nd':1.0,'Pr':1.0,'Ir':1.0,'Dy':1.0,'Er':1.0,
    'Ho':1.0,
    'R1':2.55, 'R2':2.55, 'R3':2.55, 'R4':2.55,
    'R5':2.55, 'R6':2.55, 'R7':2.55, 'R8':2.55,
    'R9':2.55, 'R10':2.55, 'R11':2.55, 'R12':2.55,
    'R13':2.55, 'R14':2.55, 'R15':2.55, 'R16':2.55,    
    'R17':2.55, 'R18':2.55, 'R19':2.55, 'R20':2.55,
    'R21':2.55, 'R22':2.55, 'R23':2.55, 'R24':2.55,    
    '?':999.0,
    }


# electronegativity difference between hydrogen and
# an h-bond donor atom.
H_DONOR_DIFFERENCE = 1.0
H_ACCEPTOR_DIFFERENCE = 1.0

# according to wikipedia, polar bonds have en-differences
# between 0.5 and 1.7.
#
# note, that with the Pauling scale, the en-difference in C-N bonds
# is 0.49, and therefore right below the lower boundary.
#
# for this reason, POLAR_BOND_MIN is set to 0.49
#
POLAR_BOND_MIN = 0.49
POLAR_BOND_MAX = 1.7

# number of outer electrons for each element.
# the Rxx elements are carbon residues as they occur in KEGG
OUTER_ELECTRONS = {
    'C':4, 'O':6, 'N':5, 'H':1,
    'R':4, 'S':6, 'P':5, 'Na':1,
    'F':7, 'Cl':7, 'Br':7, 'I':7,
    'Cu':2, 'Fe':2, 'Mg':2, 'As':9,
    'Ca':2, 'K':1, 'Li':1, 'Se':9,
    'Mn':9, 'Zn':2, 'Mo':4, 'Co':9, 
    'Ni':9, 'Cr':9, 'V':9, 'Hg':9,
    'W':9, 'Cd':9, 'Te':9, 'Pb':9,
    'Ag':1, 'Au':9, 'Pt':9, 'Si':6,
    'Sb':9, 'Bi':9, 'Ga':9, 'Sn':9,
    'Tc':9, 'B':9, 'Gd':9, 'Al':9,
    'Ba':2,'Cs':1,'Tl':3,'Pd':2,
    'Be':1,'In':1,'Ru':1,'Rh':1,'Eu':1,'La':1,
    'Ce':1,'Nd':1,'Pr':1,'Ir':1,'Dy':1,'Er':1,
    'Ho':1.0,
    'R1':4, 'R2':4, 'R3':4, 'R4':4,
    'R5':4, 'R6':4, 'R7':4, 'R8':4,
    'R9':4, 'R10':4, 'R11':4, 'R12':4,
    'R13':4, 'R14':4, 'R15':4, 'R16':4,
    'R17':4, 'R18':4, 'R19':4, 'R20':4,
    'R21':4, 'R22':4, 'R23':4, 'R24':4,
    '?':4,
    }

# sequences of bonds that are used to annotate ring structures
MAX_RING_SIZE = 7
RING_PATTERNS = {
    'C-C-C-':'cyclobutan',
    'C-C-O-C-C-':'furan',
    'C-C-C-O-C-C-':'pyran',
    'C-C-C-C-C-':'cyclopentan',
    'C-C-C-C-C-C-':'cyclohexan',    
    'C-N=C-C=C-':'pyrrol',
    'C-N-C-C=C-':'pyrrol',
    'C-N-C=C-C=':'pyrrol',
    'N-C=C-C=N-C=':'pyrimidine',
    'C-N-C=C-N=':'5-ring of purine',
    'C=N-C=C-C=C-':'pteridine',
    'C-N-C=C-C-C=':'positively_charged_pteridine',
    'N-C=C-C-N-C=':'putative_pyrimidine',
    'C-N-C-N=C-N=C-C-N=':'putative_purine',
    'C-C=C-C-C=C-C=C-C=C-C=C-C=C-C=C-C-C=C-C=':'putative_tetrapyrrol',
    'N=C-C-C-N-C-':'putative_pyrimidine_ring_of_flavin',
    'C-N-C-C=N-C-C=C-C=C-':'putative_two_rings_of_flavin',
    'C-N-C=N-C-N-C-C=N-C-C=C-C=C-':'putative_flavin',
    'C-N-C-C=N-C=':'putative_central_flavin_ring',
    'C=C-C=C-C=C-':'benzene',
    'C=N-C=C-S-':'thiazol',
    }

    
#
# atom typing scheme taken from literature
#
# The assignment is according to the table in Wang et al. 2004
# 'Development and testing of a general Amber force field'
#
#
# things not implemented:
#   biphenyl bridging carbon
#   three-member-rings
#   four-member-rings
#   inner parts of conjugated systems
#
#
# The decision tree works like
#     ("condition",
#      ("sub-condition",'result'),
#      'sub-else-result'
#      ),
#     'else-result'
#     )
#
TYPE_DECISION_TREE = (
    ("elem=='C'",(
     ("'part_of_ring' in attributes",(
      ("'sp2-hybrid' in attributes",(
       ("'aromatic_ring' in attributes",(
        ("'biphenyl_bridge_carbon' in attributes",'cp(cq)'),
        'ca')
        ),
       ("'three_member_ring' in attributes",'cu'),
       ("'four_member_ring' in attributes",'cv'),
       ("'conjugated_system' in attributes",'cc(cd)'),
       'unknown sp2-ring-carbon')
       ),
      ("'sp3-hybrid' in attributes",(
       ("'three_member_ring' in attributes",'cx'),
       ("'four_member_ring' in attributes",'cy'),
       'c3')
       ),          
      'unknown ring carbon')
      ),
     ("'sp-hybrid' in attributes",'c1'),
     ("'sp2-hybrid' in attributes",(
      ("'=O' in neighbors or '=S' in neighbors",'c'),
      ("'conjugated_system' in attributes",'ce(cf)'),
      'c2')
      ),
     ("'sp3-hybrid' in attributes",'c3'),
     'unknown carbon')
     ),
    
    ("elem=='O'",(
     ("'sp2-hybrid' in attributes",'o'),
     ("'sp3-hybrid' in attributes",(
      ("'-H' in neighbors",'o3'),
      'os')
      ),
     'unknown oxygen')
     ),

    ("elem=='N'",(
     ("'sp-hybrid' in attributes",'n1'),
     ("len(neighbors)==4",'n4'),
     ("'nitro_group' in attributes",'no'),
     ("'aromatic_ring' in attributes",'nb'),
     ("'sp2-hybrid' in attributes",(
      ("'conjugated_system' in attributes",(
       ("len(neighbors)==2",(
        ("'part_of_ring' in attributes",'nc(nd)'),
        'ne(nf)')
        ),       
       'n')
       ),
      ("len(neighbors)==3",'na'),
      'n2')
      ),    
     ("'sp3-hybrid' in attributes",(
      ("len(neighbors)==3",(
       ("atom.nb_has_feature('aromatic_ring')",'nh'),
       ("atom.nb_has_feature('amido_group')",'n'),
       'n3')
       ),
      'unknown sp3-nitrogen')
      ),
     'unknown nitrogen')
     ),
    
    ("elem=='H'",(
     ("bondtype=='HC1'",(
      ("atom.nb_has_feature('aromatic')",'ha'),
      'hc')
      ),
     ("bondtype=='HO1'",'ho'),
     ("bondtype=='HP1'",'hp'),
     ("bondtype=='HN1'",'hn'),
     ("bondtype=='HS1'",'hs'),
     'unknown hydrogenium')
     ),

    ("elem=='P'",(
     ("'aromatic_ring' in attributes",'pb'),
     ("len(neighbors)==4",'p5'),
     ("len(neighbors)==2",(
      ("'conjugated_system' in attributes",(
       ("'part_of_ring' in attributes",'pc(pd)'),
       'pe(pf)')
       ),
      'unknown 2-neighbored phosphorus')
      ),
     ("len(neighbors)==3",(
      ("'conjugated_system' in attributes",'px'),
      ("'=O' in neighbors or '=N' in neighbors or '=C' in neighbors",'p4'),
      'p3')
      ),        
     ("len(neighbors)==1",'p2'),
     'unknown phosphorus')
     ),
    
    ("elem=='S'",(
     ("len(neighbors)==1",'s2'),
     ("len(neighbors)==2",(
      ("'-H' in neighbors",'sh'),
      'ss')
      ),
     ("len(neighbors)==3",'s4'),
     ("len(neighbors)==4",'s6'),
      'unknown sulfur')
     ),
    ("elem=='F'",'f'
     ),
    ("elem=='Cl'",'cl'
     ),
    ("elem=='Br'",'br'
     ),
    ("elem=='I'",'i'
     ),
    'unk'
    )
