#!/usr/bin/env python
#
# mol_topologies.py
#
# Pattern definitions for recognizing local topologies in molecular structures.
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

import re

#
# This file contains lots of sub-molecular topologies for metabolites.
# They are conceptually similar to SMILEs, but include
# wildcard and negation operators.
#
# The patterns defined below follow a simple grammar:
# pattern == atom ( bond pattern, bond pattern , ... )
#
# where
# atom == {element, element digit, !element}
# digit == {0..9} 
# the optional digit indicating ring structures, and
#
# element == {H, He, Li, ..}
# bond == {-, =, #, ~, !-, !=, !#, !~}
# consisting of  single, double, triple and other.
# the ~ symbol references any bonds.
# the optional ! is a negative operator excluding the subsequent group

# a little dictionary containing bond valences.
BOND_SYMBOLS = {'-':1,'=':2,'#':3,'~':-1}

# a dictionary with abbreviations that are substituted, to make the patterns
# less chaotic
    
PHOSPHATE_PATTERNS = {
    'default': '~O(!~P(~O(~P)))', 
    'mono_phosphate': '-O(~P(~O(!~P),~O(!~P),~O(!~P)))', 
    'methyl_phosphonate': '-O(~P(~O(!~P),~C))', 
    'di_phosphate': '-O(~P(~O(!~P),~O(!~P),~O(~P(~O(!~P,!~C),~O(!~P,!~C),~O(!~P,!~C)))))', 
    'tri_phosphate': '-O(~P(~O(!~P),~O(!~P),~O(~P(~O(!~P),~O(!~P),~O(~P(~O(!~P),~O(!~P),~O(!~P)))))))', 
    'mod_phosphate': '~O(~P(~O(~C),~O(!~P),~O(!~P)))', 
    'mod_di_phosphate': '-O(~P(~O(!~P),~O(!~P),~O(~P(~O(!~P,~C),~O(!~P),~O(!~P)))))', 
    'amino': '~N(~C)', 
    }

SUGAR_PATTERNS = {
    'ribose':'~C9(-C(O2GROUP,-C(~O(!~P,!~V,!~C),-C(-C(PHOSPHATE),~O(~C9)))))',
    'desoxyribose':'~C9(-C(!~F,!~O,!~N,-C(-O(!~P),-C(-C(PHOSPHATE),~O(~C9)))))',
    'alpha-desoxyribose':'~C9(-C(!~O,!~N,-C(-O(!~P),-C(-C(PHOSPHATE),~O(~C9)))))',
    'tetrose':'~C9(-C(O2GROUP,-C(PHOSPHATE,-C(!-C,-O(-C9)))))',
    'arabinose':'~C9(-C(O2GROUP,-C(~O(!~P,!~V,!~C),-C(-C(PHOSPHATE),~O(~C9)))))',
    'ribose_3p':'~C9(-C(O2GROUP,-C(~O(~P),-C(-C(PHOSPHATE),~O(~C9)))))',
    'desoxyribose_3p':'~C9(-C(!~O,-C(~O(~P),-C(-C(PHOSPHATE),~O(~C9)))))',
    'ribose_23p':'~C9(-C(~O(~P(~O8)),-C(~O8(~P),-C(-C(PHOSPHATE),~O(~C9)))))',
    'ribose_3meo':'~C9(-C(O2GROUP,-C(~O(~C(-C(-O(-C)))),-C(-C(PHOSPHATE),~O(~C9)))))',
    'ribose_23vanadate':'~C9(-C(~O(~V(~O8)),-C(~O8,-C(-C(PHOSPHATE),~O(~C9)))))',
    'phosphono_ribose':'~C9(-C(~O(~C(~O8,~P(~O,~O))),-C(~O8,-C(-C(PHOSPHATE),~O(~C9)))))',
    }

O2GROUP_PATTERNS = {
    'oh': '~O(!~C)', 
    'methyl-o': '~O(~C(!~C,!~O))', 
    'mcarb':'~O(~C(~O,~N(~C)))', 
    'fluoride':'~F', 
    '2n':'~N(!~C)', 
    '2n-acetyl':'~N(~C(~O,~C))', 
    'ethyl':'~O(~C(~C))', 
    }

BASE_PATTERNS = {
    'adenine':'N1(~C(!-C,!~N,!~O,~N(~C(~N(!~C),~C2(~C(~N1,~N(SUGAR,~C(~N(~C2)))))),!-C,!~N)))',
    'uracil': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(=C(~N1,!~C,!~O),!~C,!~O,!~F)),!~C)))',
    'cytosine':'N1(SUGAR,~C(~O,~N(!-C,~C(~N(!-C),~C(~C(~N1),!~C,!~N,!~F)))))',
    'guanine': 'N1(~C(~N(~C(~O(!~C),~C2(~C(~N1,~N(SUGAR,~C(!~C,!~O,!~N(~C,~C),~N(~C2,!-C)))))),!-C),~N(!~C)))',
    'thymine': 'N1(SUGAR,~C(~O,~N(!~C,~C(~O,~C(-C(!~O,!~N,!~C),~C(~N1))))))',        
    
    # natural modifications
    'inosine': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(~C2)))))),!-C),!~N))',
    'pseudouracil':'C1(SUGAR,~C(~O,~N(!-C,~C(~O,~N(~C(~C1),!-C)))))',
    'isoguanine': 'N1(~C(~N(~C(~N(!~C),~C2(~C(~N1,~N(SUGAR,~C(!~C,!~O,!~N(~C,~C),~N(~C2,!-C)))))),!-C),~O(!~C)))',

    'Arp':'N1(~C(~N(~C(~N,~C2(~C(~N1,~N(~C(~N(~C2)),-C9(-C(-O(-C8(-C(-O(!-C),-C(-O,-C(-C(~O(~P)),-O(-C8)))))),-C(-O,-C(-O(-C9))))))))))))',
    'D': 'N1(~C(~O,~N(~C(~O,-C(-C(-N1,!~C),!-C,!~N)))),-C,-C)',
    'Grp': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(-C3(-C(-O(-C4(-C(-O,-C(-O,-C(~C(~O(~P)),-O(-C4)))))),-C(-O,-C(-O(-C3))))),~C(~N(~C2))))))),~N))',
    'P': 'N1(~C3(~N(~C(~C(~C(~N(~C3)))),~C(~O,~C2(~C(~N1,~N(~C(~N(~C2,!-C)))))))))',
    'QtRNA': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(-C(-N(-C3(-C(!-O,~C(!-O,~C(-O(!-C),-C(-C3,-O(!-C)))))))),~C2))))))),~N))',
    'R': 'N1(~C(~N(~C(-N(~C(~C(~O),~C(~C(~C(~C))))),~C2(~C(~N1,~N(~C(~N(~C2)))))),!-C)))',
    
    'ac6A':'N1(~C(~N(~C(~N(-C(~O,-C)),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'ac4C':'N1(~C(~O,~N(~C(~N(-C(~O,-C)),~C(~C(~N1))))),SUGAR)',
    'acp3U': 'N1(~C(~O,~N(-C(-C(-C(-N,-C(~O,~O)))),~C(~O,~C(~C(~N1))))),SUGAR)',
    'G+': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(~C2,-C(~N,~N)))))))),~N))',
    'chm5U':'N1(~C(~O,~N(~C(~O,~C(-C(-O,-C(~O(!~C),~O(!~C))),~C(~N1))))))',
    'cm5U': 'N1(~C(~O,~N(~C(~O,~C(-C(-C(~O(!~C),~O(!~C)),!-O),~C(~N1))))))',
    'cm5s2U': 'N1(~C(~S,~N(~C(~O,~C(-C(-C(~O(!~C),~O(!~C))),~C(~N1))))))',
    'cmnm5U': 'N1(~C(~O,~N(~C(~O,~C(-C(-N(-C(-C(~O,~O)))),~C(~N1))))),SUGAR)', 
    'cmnm5s2U':'N1(~C(~S,~N(~C(~O,~C(-C(-N(-C(-C(~O,~O)))),~C(~N1))))))',
    'cmnm5se2U':'N1(~C(~Se,~N(~C(~O,~C(-C(-N(-C(-C(~O,~O)))),~C(~N1))))))',
    'cmo5U':'N1(~C(~O,~N(~C(~O,~C(~O(~C(-C(~O(!~C),~O(!~C)))),~C(~N1))))))',
    'f5C': 'N1(~C(~O,~N(~C(~N,~C(~C(=O),~C(~N1))))),SUGAR)',
    'f5U': 'N1(~C(~O,~N(~C(~O,~C(~C(~O(!~C),!-C),~C(~N1))))),SUGAR)',    
    'f5s2C': 'N1(~C(~S,~N(~C(~N,~C(~C(~O),~C(~N1))))))',
    'f5s2U': 'N1(~C(~S,~N(~C(~O,~C(~C(~O),~C(~N1))))))',
    'f5se2C':'N1(~C(~Se,~N(~C(~N,~C(~C(~O),~C(~N1))))))',    
    'f5se2U':'N1(~C(~Se,~N(~C(~O,~C(~C(~O),~C(~N1))))))',
    'g6A':'N1(~C(~N(~C(~N(!-C,~C(~O,~N(-C(-C(~O,~O),!-C)))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'galQtRNA':'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(-C(-N(-C3(-C(~C(-C(-O,-C(-C3,-O(-C4(-C(-O,-C(-O,-C(-O,-C(-C(-O),-O(-C4)))))))))))))),~C2))))))),~N))',
    'gluQtRNA1':'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(-C(-N(-C3(-C(~C(-C(~O(~C(~O,-C(-N,-C(-C(-C(~O,~O)))))),-C(-C3,-O))))))),~C2))))))),~N))',
    'gluQtRNA2':'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(-C(-N(-C3(-C(~C(-C(-O,-C(-C3,~O(~C(~O,-C(-N,-C(-C(-C(~O,~O))))))))))))),~C2))))))),~N))',
    'hm5C': 'N1(~C(~O,~N(~C(~N,~C(-C(-O),~C(~N1))))))',
    'hn6A': 'N1(~C(!-S,~N(~C(~N(-C(~O,-N(-C(-C(-O,-C(-C)),-C(~O,~O))))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'ho5U': 'N1(~C(~O,~N(~C(~O,~C(~O(!~C),~C(~N1))))))',
    'i6A': 'N1(~C(~N(~C(~N(-C(-C(~C(-C(!~O),-C(!~O))))),~C2(~C(~N1,~N(~C(~N(~C2))))))),!-S))',
    'imG2': 'N1(~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(-C,~C(-C(!-C),~N3)))),~C,!-C)',
    'imG': 'N1(-C(!~C),~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(-C,~C(!-C,~N3)))))',
    'imG-14': 'N1(~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(-C,~C(!-C,~N3)))),~C,!-C)',
    'mimG': 'N1(-C(!~C),~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(-C,~C(-C(!-C),~N3)))))',
    'io6A': 'N1(~C(!-S,~N(~C(~N(-C(-C(~C(-C(-O),-C)))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'k2C': 'N1(~C(~N(-C(-C(-C(-C(-C(-C(~O,~O),-N)))))),~N(~C(~N,~C(~C(~N1,!~N),!~N)))))',    
    'm1A': 'N1(~C(~N(-C,~C(~N(!-C),~C2(~C(~N1,~N(SUGAR,~C(~N(~C2)))))))))',
    'm1G': 'N1(~C(~N(-C(!~C),~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(~C2))))))),~N))',
    'm1I': 'N1(~C(~N(-C(!~C),~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(~C2))))))),!~N))',    
    'm1Y': 'C1(~C(~O,~N(!-C,~C(~O,~N(-C(!-C,!-O),~C(~C1))))))',
    'm1acp3Y': 'C1(SUGAR,~C(~O,~N(-C(-C(-C(-C(~O,~O),-N))),~C(~O,~N(-C,~C(~C1))))))',    
    'm22G': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(!-C,~C2))))))),~N(-C(!~C),-C)))',    
    'm2A': 'N1(~C(~N(~C(~N,~C2(~C(~N1,~N(~C(~N(~C2))))))),-C))',
    'm2G': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(!-C,~C2))))))),~N(-C(!~C),!-C)))',    
    'm3C': 'N1(SUGAR,~C(~O,~N(-C,~C(~N,~C(~C(~N1))))))',
    'm3U': 'N1(SUGAR,~C(~O,~N(-C(!-C),~C(~O,~C(~C(~N1))))))',
    'm3Y': 'C1(SUGAR,~C(~O,~N(-C(!-C),~C(~O,~N(~C(~C1))))))',
    'm44C': 'N1(SUGAR,~C(~O,~N(~C(~N(-C,-C),~C(~C(~N1))))))',
    'm4C': 'N1(SUGAR,~C(~O,~N(~C(~N(-C(!~O,!-C),!-C),~C(~C(~N1))))))',
    'm5C': 'N1(~C(~O,~N(~C(~N,~C(~C(!~O),~C(~N1))))),SUGAR)',    
    'm5D': 'N1(~C(~O,~N(~C(~O,-C(-C(!~O,!~N),-C(-N1))))),-C,-C))',
    'm5s2U': 'N1(~C(~S,~N(~C(~O,~C(-C(!-C,!~N,!~O),~C(~N1))))))',
    'm66A': 'N1(~C(~N(~C(~N(-C(!~C,!~N,!~O),-C(!~C,!~N,!~O)),~C2(~C(~N1,~N(SUGAR,~C(~N(~C2)))))))))',
    'm6A': 'N1(~C(!-S,~N(~C(~N(-C(!~C,!~N),!~C),~C2(~C(~N1,~N(SUGAR,~C(~N(~C2)))))))))',
    'm6t6A': 'N1(~C(!-S,~N(~C(~N(-C,~C(~O,~N(-C(-C(~O,~O),-C(-C,-O))))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'm7G': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(-C(!~C),~C2))))))),~N(!-C)))',
    'm7m22G': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~N(-C,~C2))))))),~N(-C,-C)))', 
    'm7m2G': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~N(-C,~C2))))))),~N(-C,!-C)))',
    'manQtRNA': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(-C(-N(-C3(-C(~C(-C(-O,-C(-C3,-O(-C4(-C(-O,-C(-O,-C(-O,-C(-C(-O),-O(-C4)))))))))))))),~C2))))))),~N))',
    'mchm5U': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(-C(-O,-C(~O,~O(~C))),~C(~N1))))))',
    'mcm5U'  : 'N1(SUGAR,~C(~O,~N(~C(~O,~C(-C(!~O,-C(~O,~O(~C))),~C(~N1))))))',
    'mcm5s2U': 'N1(~C(~S,~N(~C(~O,~C(-C(-C(~O,~O(~C))),~C(~N1))))))',
    'mcmo5U': 'N1(~C(~O,~N(~C(~O,~C(~O(~C(-C(~O,~O(~C)))),~C(~N1))))))',
    'mnm5U': 'N1(~C(~O,~N(~C(~O,~C(-C(-N(-C(!~C))),~C(~N1))))))',
    'mnm5s2U': 'N1(~C(~S,~N(~C(~O,~C(-C(-N(-C(!-C))),~C(~N1))))))',
    'mnm5se2U': 'N1(~C(~Se,~N(~C(~O,~C(-C(-N(-C(!-C))),~C(~N1))))))',
    'mo5U': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(~O(-C(!~C)),~C(~N1))))))',
    'mh5U': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(~C(~O(!~C)),~C(~N1))))))',
    'ms2hn6A': 'N1(~C(-S(-C),~N(~C(~N(-C(~O,-N(-C(-C(-O,-C(-C)),-C(~O,~O))))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'ms2i6A': 'N1(~C(-S(-C),~N(~C(~N(-C(-C(~C(-C(!-O),-C(!-O))))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'ms2io6A': 'N1(~C(-S(-C),~N(~C(~N(-C(-C(~C(-C(-O),-C)))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'ms2m6A': 'N1(~C(-S(-C),~N(~C(~N(-C(!~C,!-N),!-C),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'ms2t6A': 'N1(~C(-S(-C),~N(~C(~N(-C(~O,-N(-C(-C(~O,~O),-C(-C(!-C),-O))))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',    
    'ncm5U': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(-C(-C(~O,-N)),~C(~N1))))))',
    'nm5U': 'N1(~C(~O,~N(~C(~O,~C(-C(-N(!-C)),~C(~N1))))))',
    'nm5s2U': 'N1(~C(~S,~N(~C(~O,~C(-C(-N(!-C)),~C(~N1))))))',
    'nm5se2U': 'N1(~C(~Se,~N(~C(~O,~C(-C(-N(!-C)),~C(~N1))))))',    
    'o2yW': 'N1(~C,~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(~C,~C(-C(-C(-C(-C(~O(~C),~O),~N(~C(~O(~C),~O))),~O(~O))),~N3)))))',
    'oHyW': 'N1(~C,~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(~C,~C(-C(-C(-C(-C(~O(~C),~O),~N(~C(~O(~C),~O))),~O(!~O))),~N3)))))',
    "oHyWx": 'N1(-C,~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(~C,~C(-C(-C(-C(-C(~O(~C),~O),~N(!~C)),~O)),~N3)))))',
    'oQtRNA': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(-C(-N(-C3(-C4(~C(~O(~C4),-C(-O,-C(-C3,-O))))))),~C2))))))),~N))',
    'preQ0tRNA': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(~C(#N(!~C,!~O),!~C),~C2))))))),~N))',
    'preQ1tRNA': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(~C(~C(-C(-N(!-C),!~N),~C2))))))),~N))',    
    's2C': 'N1(SUGAR,~C(~S(!-C),~N(~C(~N,~C(~C(~N1,!~C),!~C)))))',
    's2U': 'N1(SUGAR,~C(~S,~N(~C(~O,~C(!~C,~C(~N1))))))',
    's4U': 'N1(~C(~O,~N(~C(~S,~C(~C(~N1))))))',
    'se2U': 'N1(~C(~Se,~N(~C(~O,~C(!~C,~C(~N1))))))', 
    't6A':'N1(~C(!-S,~N(~C(~N(!-C,~C(~O,~N(~C(~C(~O,~O),~C(~C(!-C),~O))))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    'tm5U': 'N1(~C(~O,~N(~C(~O,~C(-C(-N(-C(-C(~S(~O,~O,~O))))),~C(~N1))))))',
    'tm5s2U': 'N1(~C(~S,~N(~C(~O,~C(-C(-N(-C(-C(~S(~O,~O,~O))))),~C(~N1))))))',
    'yW': 'N1(~C(!~C),~C(~N3(~C(~O,~C2(~C(~N1,~N(~C(~N(~C2))))))),~N(~C(-C,~C(-C(-C(!-O,-C(~N(~C(~O,~O(~C))),-C(~O,~O(~C))))),~N3)))))',

    # other base modifications found
    'ms2i6Aiso': 'N1(~C(-S(-C),~N(~C(~N(-C(=C(-C(-C(!-O),-C)))),~C2(~C(~N1,~N(~C(~N(~C2)))))))))',
    '8fG': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(~C,~C),~N(~C2,!-C)))))),!-C),~N(!~C)))',
    '8oG': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~O(!~C),~N(~C2,!-C)))))),!-C),~N(!~C)))',
    'm6G': 'N1(~C(~N(~C(~O(-C),~C2(~C(~N1,~N(SUGAR,~C(~N(~C2))))))),~N(!-C)))',
    'm8G': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~C,~N(~C2,!-C)))))),!-C),~N(!~C)))',
    'N72G': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(!-C,~C(!~O,!~N(~C,~C),~N(~C2,-C)))))),!-C),~N(!~C)))',
    '5nitroC': 'N1(~C(~O,~N(~C(~N,~C(~N(~O,~O),~C(~N1))))),SUGAR)',    
    '5nitroU': 'N1(~C(~O,~N(~C(~O,~C(~N(~O,~O),~C(~N1))))),SUGAR)',    
    'N2A': 'N1(~C(~N(~C(~N,~C2(~C(~N1,~N(SUGAR,~C(~N(~C2))))))),~N))',
    '5propU': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(~C(~C(~C(!~C,!~N))),~C(~N1))))))',        
    '5propnU': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(~C(~C(~C(!~C,~N))),~C(~N1))))))',        
    'ethC':'N1(SUGAR,~C(~O,~N2(~C(~N(~C(~C(~N2))),~C(~C(~N1),!~C,!~N)))))',
    '5fluoroU': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(~F,~C(~N1))))))',        
    '5fluoroC': 'N1(SUGAR,~C(~O,~N(~C(~N,~C(~F,~C(~N1))))))',        
    'dhOro': 'N1(~C(~O,~N(~C(~O,-C(-C(-N1,-C(~O,~O)),!-C,!~N)))),-C,-C)',
    '7fluorobenzylG': 'N1(~C(~N(~C(~O,~C2(~C(~N1,~N(SUGAR,~C(~N(~C(~C3(~C(~C(~C(~C(~C(~C3))))))),~C2))))))),~N(!-C)))',
    'e3m5U': 'N1(SUGAR,~C(~O,~N(~C(~C(!~C)),~C(~O,~C(-C(!~O,!~N,!~C),~C(~N1))))))',        
    'o6U': 'N1(SUGAR,~C(~O,~N(~C(~O,~C(=C(~N1,!~C,~O),!~C,!~O,!~F)),!~C)))',
    }

OLD_MODS = (
    ('purine','N1(~C(~N(~C(~C2(~C(~N1,~N(~C(~N(~C2)))))))))'),
    ('pyrimidine','N1(~C(~N(~C(~C(~C(~N1,!~N)),!~N))))'),
    ('riboside','C1(-C(-C(-C(-O(-C1)),-O),-O),-N)'),
    ('desoxyriboside','C1(-C(-C(-C(-O(-C1)),-O),!-O),-N)'),
    ('phosphate','P(~O,~O,-O(-C))'),

)




PATTERN_ABBREVIATIONS = (
    ('cbx','C(=O,-O)'),
    ('meth','C(-H,-H,-H)'),
    )

class Pattern(object):
    """
    Pre-parsed topology pattern.
    """
    def __init__(self, name, pattern):
        self.name = name
        for abbrev, full in PATTERN_ABBREVIATIONS:
            pattern = re.sub(abbrev, full,pattern)
        self.pattern = pattern
        self.topology = self._create_topology(self.pattern)
        self.sumdict = self._create_sumdict(self.topology)

    def _create_topology(self, pattern, bool=1,valence=0):
        """
        Parses a topology string pattern into something that the
        recursive matching algorithm in topology_matcher.match_node() can
        handle quickly.

        Returns a NODE consisting of (ELEM,SUBNODES,BOOL,ID,VALENCE), where
        ELEM is a string containing the element at this node,
        SUBNODES is a list of NODE tuples
        (describing bond order and bound atom),
        BOOL is either 0 or 1 and tells whether the topology behind
        this node is required (1) or excluded (0).
        ID is an index number used to match atoms that occur in an
        expression more than once. It is normally zero.
        VALENCE is an integer describing how this atom should be bonded
        by its predecessor
           1-single,2-double,3-triple,-1-arbitrary,0-toplevel node.

        The parsing algorithm also works recursively.
        Have fun debugging.
        """
        jmax = len(pattern)
            
        elem = ""
        j = 0
        index = ""
        while j<jmax and pattern[j] not in ('(',')',','):
            # handle digit indices
            if pattern[j] in '0123456789': index += pattern[j]
            else: elem += pattern[j]
            j += 1
        # all digits following elements are considered ID's
        if index: index = int(index)
        else: index = 0
        
        # string end, ',' and ')' end this node, jump back.
        if j==jmax or pattern[j] in (',',')'):
            subnodes = []
        else:
            # upon '(', sub-nodes need to be parsed, until the parentheses close.
            j += 1
            subnodes = []        
            # collect all comma separated nodes until parentheses close
            node_string = ""
            parenthese_level = 0
            while j < jmax and parenthese_level >= 0:
                char = pattern[j]
                if char == '(': parenthese_level += 1
                if char == ')': parenthese_level -= 1
                if (char == ',' and parenthese_level == 0) or (char == ')' and parenthese_level == -1):
                    # look for negation '!' symbol at string start.
                    # '!' excludes the entire following expression
                    if node_string[0] == '!': boolean = 0
                    else: boolean = 1
                    # next character must be the valence symbol
                    next_valence = BOND_SYMBOLS[node_string[1-boolean]]

                    # everything else is the sub-pattern
                    # that needs to be parsed by the recursive procedure call.
                    node_string = node_string[2-boolean:]
                    subnodes.append(self._create_topology(node_string,boolean,next_valence))
                    node_string = ""
                else:
                    node_string += pattern[j]
                j += 1
                
            if parenthese_level != -1:
                print(("MISMATCHED PARENTHESES IN PATTERN:",pattern))

        return (elem,subnodes,bool,index,valence)

    def get_sumdict_node(self, sumdict,node,circles):
        """Recursive"""
        elem,subnodes,boole,idd,valence = node
        if boole==0: return
        if idd!=0:
            if idd in circles: return
            else:
                circles[idd] = True
        sumdict.setdefault(elem,0)
        sumdict[elem] += 1
        for sn in subnodes:
            self.get_sumdict_node(sumdict,sn,circles)
                
    def _create_sumdict(self, topology):
        sumdict = {}
        self.get_sumdict_node(sumdict,topology,{})
        return sumdict

    def __repr__(self):
        return self.pattern


class NucleotidePattern(Pattern):
    def __init__(self, name, base='adenine', sugar='ribose', phosphate='not_specified', o2group='oh'):
        pattern = BASE_PATTERNS[base]
        pattern = pattern.replace('SUGAR', SUGAR_PATTERNS[sugar])
        pattern = pattern.replace('PHOSPHATE', PHOSPHATE_PATTERNS[phosphate])
        pattern = pattern.replace('O2GROUP', O2GROUP_PATTERNS[o2group])
        pattern = pattern.replace('()', '')
        Pattern.__init__(self, name, pattern)
    

def read_nucleotide_topologies(filename):
    """Prepares a list of nucleotide patterns."""
    result = []
    for line in open(filename):
        tokens = line.strip().split()
        if not line.startswith('#') and len(tokens)>1:
            name, base, sugar, phosphate, o2group = tokens
            pattern = NucleotidePattern(name, base=base, sugar=sugar, \
                                        phosphate=phosphate, o2group=o2group)
            result.append(pattern)
    return result
    


# this list consists of all organic chemistry groups i could find,
# that could play a role in cellular metabolisms
FUNCTIONAL_GROUP_PATTERNS = (
    ('h_hydrogen_h2','H(-H)'),
    ('h_water','H(-O(-H))'),
    ('h_other_hydrogen','H(!-H,!-O(-H))'),
    ('o_oxygen','O(=O)'),
    ('o_peroxide','O(-O)'),
    ('o_water','O(-H,-H)'),
    ('o_carbon_dioxide','O(=C(=O))'),
    ('c_carbon_dioxide','C(=O,=O)'),
    ('n_nitrate','N(=O,=O,-O(-H))'),
    ('n_nitrite','N(=O,!=O,-O(-H))'),
    ('chinone','C1(=O,-C(=C(-C(=O,-C(=C(-C1))))))'),
    ('hydrochinone','C1(-O(-H),-C(=C(-C(-O(-H),=C(-C(=C1))))))'),
    ('o_methoxy','O(-methyl,-C(!-O,!=O,!-N,!=N))'),
    ('c_methoxy','C(-H,-H,-H,-O(-C(!-O,!=O,!-N,!=N)))'),
    ('c_with_methoxy','C(-O(-methyl),!-O,!=O,!-N,!=N)'),
    ('n_ammonia','N(-H,-H,-H)'),
    ('n_primary_amine','N(-H,-H,-C(!-O,!=O))'),
    ('n_secondary_amine','N(-H,-C(!-O,!=O),-C(!-O,!=O))'),
    ('n_tertiary_amine','N(-C(!-O,!=O),-C(!-O,!=O),-C(!-O,!=O),!=C)'),
    ('n_quarternary_amine','N(-C(!-O,!=O),-C(!-O,!=O),-C(!-O,!=O),-C(!-O,!=O))'),
    ('c_primary_amine','C(-N(-H,-H,!-C,!=C),!-O,!=O)'),
    ('c_secondary_amine','C(-N(-H,-C(!-O,!=O),!-C,!=C),!-O,!=O)'),
    ('c_tertiary_amine','C(-N(-C(!-O,!=O),-C(!-O,!=O),!-C),!-O,!=O)'),
    ('c_quarternary_amine','C(-N(-C(!-O,!=O),-C(!-O,!=O),-C(!-O,!=O)),!-O,!=O)'),
    ('c_methyl','C(-H,-H,-H)'),
    ('c_secondary','C(-H,-H,!-H,!~O,!~N,!~S)'),
    ('c_tertiary','C(-H,!-H,!~O,!~N,!~S)'),
    ('c_quarternary','C(!-H,!~O,!~N,!~S)'),
    ('c_ethyl_c1','C(-H,-H,-methyl)'),
    ('c_ethyl_c2','C(-H,-H,-H,-C(-H,-H))'),
    ('c_propyl_c1','C(-H,-H,-C(-H,-H,-methyl))'),
    ('c_propyl_c2','C(-H,-H,-methyl,-C(-H,-H)))'),
    ('c_propyl_c3','C(-H,-H,-H,-C(-H,-H),-C(-H,-H)))'),
    ('c_butyl_c1','C(-H,-H,-C(-H,-H,-C(-H,-H,-methyl)))'),
    ('c_butyl_c2','C(-H,-H,-C(-H,-H,-methyl),-C(-H,-H))'),
    ('c_butyl_c3','C(-H,-H,-methyl,-C(-H,-H,-C(-H,-H)))'),
    ('c_butyl_c4','C(-H,-H,-H,-C(-H,-H,-C(-H,-H,-C(-H,-H))))'),
    ('c_tert_butyl','C(-methyl,-methyl,-methyl)'),
    ('c_isopropyl','C(-methyl,-methyl,-H)'),
    ('c_methylene','C(-H,-H,!-H)'),
    ('c_cyanid','C(#N)'),
    ('n_cyanid','N(#C)'),
    ('o_nitro','O(-N(=O,-C),!-C,!-N,!-H)'),
    ('o_nitro','O(=N(-O,-C),!-C,!-N,!-H)'),    
    ('n_nitro','N(-O(!-C,!-N,!-H),=O,-C)'),
    ('c_nitro','C(-N(-O(!-C,!-N,!-H),=O))'),
    ('c_amid','C(=O,-N(!=C))'),
    ('o_amid','O(=C(-N(!=C)))'),
    ('n_amid','N(-C(=O)!=C))'),
    ('c_thioether','C(-S(-C(!=O)),!=O)'),
    ('s_thioether','S(-C(!=O),-C(!=O))'),
    ('c_thioester','C(=O,-S(-C))'),
    ('s_thioester','S(-C(=O),-C)'),
    ('o_thioester','O(=C(-S(-C)))'),    
    ('c_primary_thiol','C(-S(-H),-H,-H)'),
    ('c_secondary_thiol','C(-S(-H),-H,!-H)'),
    ('c_tertiary_thiol','C(-S(-H),!-H,!-H)'),
    ('s_primary_thiol','S(-H,-C(-H,-H))'),
    ('s_secondary_thiol','S(-H,-C(-H,!-H))'),
    ('s_tertiary_thiol','S(-H,-C(!-H,!-H))'),
    ('c_carbonyl','C(=O,!-O(-H),!-S,!-N,!=O)'),
    ('c_aldehyde','C(=O,-H,!-N,!-S)'),
    ('c_ketone','C(=O,!-O,!=O,!-N,!=N,!-H)'),
    ('o_carbonyl','O(=C(!-O,!-N,!-S,!=O,!=N))'),
    ('c_ester','C(=O,-C,-O(-C(!=O,!-O)))'),
    ('c_half_acetal','C(-H,-O(-C(!=O)),-O(-H))'),    
    ('c_half_ketal','C(!-H,-O(-C(!=O)),-O(-H))'),
    ('c_full_acetal','C(-O(-C(!=O)),-O(-C(!=O)))'),
    ('o_acetal','O(-C(!=O),-C(!=O,-O))'),
    ('c_primary_alcohol','C(-O(-H),-H,-H)'),
    ('c_secondary_alcohol','C(-O(-H),-H,-C,-C)'),
    ('c_tertiary_alcohol','C(-O(-H),-C,-C,-C)'),
    ('o_hydroxyl','O(-H,-C(!=O,!-O,!-N,!=N,!-S))'),
    ('o_ether','O(-C(!-O,!=O,!-N,!=N,!-S),-C(!-O,!=O,!-N,!=N,!-S))'),
    ('o_carboxyl','O(-H,-C(=O,!-O))'),
    ('o_carboxyl','O(=C(-O(-H),!-O))'),
    ('c_carboxyl','C(=O,-O(-H),!-O))'),
    ('o_carbonic_acid_anhydride','O(-C(=O,!-O),-C(=O,!-O))'),
    ('c_carbonic_acid_anhydride','C(=O,-O(-C(=O,!-O)),!-O)'),
    ('c_carbonic_acid_chloride','C(=O,-Cl)'),
    ('o_carbonic_acid_chloride','O(=C(-Cl))'),
    ('cl_carbonic_acid_chloride','Cl(-C(=O))'),
    ('c_enol','C(=C,-O(-H))'),
    ('o_enol','O(-H,-C(=C))'),
    ('s_sulfonamide','S(=O,=O,-N(-C),-C)'),
    ('n_sulfonamide','N(-C,-S(=O,=O,-C))'),
    ('c_sulfonamide','C(-N(-S(=O,=O,-C)))'),
    ('o_sulfonamide','O(=S(=O,-C,-N(-C)))'),
    ('c_phosphoester','C(-O(-P(=O,-O,-O)),!=O,!-N)'),
    ('p_orthophosphate','P(-O(!-C,!-N,!-P),-O(!-C,!-N,!-P),-O(!-C,!-N,!-P),=O)'),    
    ('p_phosphodiester','P(=O,-O,-O(-C(!-N,!=O)),-O(-C(!-N,!=O)))'),
    ('p_phosphoester','P(=O,-O(-C(!-N,!=O)),-O,-O)'),
    ('p_phosphoanhydrid','P(=O,-O,-O,-O(-P(=O,-O,-O)))'),
    ('p_phosphodianhydrid','P(=O,-O,-O(-P(=O,-O,-O)),-O(-P(=O,-O,-O)))'),
    ('beta-carboxy-carboxylate','C(-O,=O,-C(-cbx))'),
    ('s_thiobridge','-C(-S(-S(-C)))'),
    ('hydroxy_benzene','C1(=C(-C(=C(-C(=C(-C1))))),-O(-H))'),
    ('phenolic_ester','C1(=C(-C(=C(-C(=C(-C1))))),-O(-C(=O)))'),
    ('phenolic_phosphate_ester','C1(=C(-C(=C(-C(=C(-C1))))),-O(-P(=O)))'),
    ('n_imine','N(!-O,-C(!-O,!=O),=C)'),
    ('c_imine','C(=N(=C,!-O),!-O,!=O)'),
    ('s_methylthioether','S(-C,-methyl)'),
    ('o-glycoside','C(-O(-C(-O(-C(!~O)),!-O)))'),
    ('n-glycoside','C(-O(-C(-N(-C(!~O)),!-O)))'),
    ('aminoacyl','C(-N(-H),-H,-C(=O))'),
    ('glycyl','C(-N(-H),-H,-H,-C(=O))'),    
    ('alanyl','C(-N(-H),-H,-methyl,-C(=O))'),
    ('valyl','C(-N(-H),-H,-C(-methyl,-methyl,-H),-C(=O))'),
    ('leucyl','C(-N(-H),-H,-C(-H,-H,-C(-H,-methyl,-methyl)),-C(=O))'),
    ('isolecyl','C(-N(-H),-H,-C(-methyl,-C(-H,-H,-methyl)),-C(=O))'),    
    ('cysteinyl','C(-N(-H),-H,-C(-S,-H,-H),-C(=O))'),
    ('phenylalanyl','C(-N(-H),-H,-C(-H,-H,-C1(=C(-H,-C(-H,=C(-H,-C(-H,=C(-H,-C1))))))),-C(=O))'),
    ('tyrosyl','C(-N(-H),-H,-C(-H,-H,-C1(=C(-H,-C(-H,=C(-O,-C(-H,=C(-H,-C1))))))),-C(=O))'),
    ('seryl','C(-N(-H),-H,-C(-H,-H,-O),-C(=O))'),
    ('threonyl','C(-N(-H),-H,-C(-H,-O,-methyl),-C(=O))'),    
    ('glutamyl','C(-N(-H),-H,-C(-H,-H,-C(-H,-H,-C(=O))),-C(=O))'),
    ('aspartyl','C(-N(-H),-H,-C(-H,-H,-C(=O)),-C(=O))'),
    ('prolyl','C1(-N(-C(-H,-H,-C(-H,-H,-C(-H,-H-C1)))),-H,-C(=O))'),
    ('histidyl','C(-N(-H),-H,-C(-H,-H,-C1(~C(~N(~C(~N(~C1)))))),-C(=O))'),    
    ('arginyl','C(-N(-H),-H,-C(-H,-H,-C(-H,-H,-C(-H,-H,-N(-C(=N,-N))))),-C(=O))'),    
    ('methionyl','C(-N(-H),-H,-C(-H,-H,-C(-H,-H,-S(-C(-H,-H)))),-C(=O))'),
    ('lysyl','C(-N(-H),-H,-C(-H,-H,-C(-H,-H,-C(-H,-H,-C(-H,-H,-N)))),-C(=O))'),    
    ('tryptophanyl','C(-N(-H),-H,-C(-H,-H,-C1(~C2(~C3(~C(~C(~C(~C(~C2)))),~N(~C(~C1)))))),-C(=O))'),
    ('pyranose','C7(-C(-C(-C(-C(-O(-C7)),-O),-O),-O))'),
    ('furanose','C8(-C(-O,-C(-O,-C(-O(-C8)))))'),
    ('sterol','C1(~C2(-C3(-C4(-C5(-C6(-C(-C(-C(-C5))),-C(-C(-C3)))),~C(~C(~C1)))),~C(~C(~C(~C(~C1))))))'),
    ('c_trifluoromethyl','C(-F,-F,-F)'),
    ('c_trichloromethyl','C(-Cl,-Cl,-Cl)'),
    )

# This is an extra set of topologies that are used for finding additional groups in alpha position
ALPHA_GROUP_PATTERNS = (
    ('peroxide','O(-O)'),
    ('trifluoromethyl','C(-F,-F,-F)'),
    ('trichloromethyl','C(-Cl,-Cl,-Cl)'),
    ('hydroxyl','C(-O(-H),!~O,!~N)'),
    ('carbonyl','C(=O,!~O,!~S,!~N)'),
    ('carboxy','C(-O,=O)'),
    ('amino','C(-N),!=O,!-N'),
    ('alhdehyde','C(=O,-H)'),
    ('keto','C(=O,!-H,-C,!-O)'),
    ('diol',  'C(-O(-H),!~O,-C(-O(-H),!~O))'),
    ('triol','C(-O(-H),!~O,-C(-O(-H),!~O,-C(-O(-H),!~O)))'),
    ('ester','C(=O,-O(-C(!~O)))'),
    ('reverse_ester','C(!=O,-O(-C(=O)))'),    
    ('half_acetal','C(-H,-O(-C(!=O)),-O(-H))'),    
    ('half_ketal','C(!-H,-O(-C(!=O)),-O(-H))'),
    ('full_acetal','C(-O(-C(!=O)),-O(-C(!=O)))'),
    ('primary_alcohol','C(-O(-H),-H,-H)'),
    ('secondary_alcohol','C(-O(-H),-H,-C,-C)'),
    ('tertiary_alcohol','C(-O(-H),-C,-C,-C)'),
    ('ether','C(-O(-C(!~O,!~N,!~S)),!~O,!~N,!~S))'),
    ('carbonic_acid_anhydride','C(=O,-O(-C(=O,!-O)),!-O)'),
    ('carbonic_acid_chloride','C(=O,-Cl)'),
    ('enol','C(=C,-O(-H))'),
    ('reverse_enol','C(=C(-O(-H)))'),
    ('en','C(=C)'),
    ('conjugated','C(=C(-C(=C)))'),
    ('primary_amine','C(-N(-H,-H,!-C,!=C),!-O,!=O)'),
    ('secondary_amine','C(-N(-H,-C(!-O,!=O),!-C,!=C),!-O,!=O)'),
    ('tertiary_amine','C(-N(-C(!-O,!=O),-C(!-O,!=O),!-C),!-O,!=O)'),
    ('quarternary_amine','C(-N(-C(!-O,!=O),-C(!-O,!=O),-C(!-O,!=O)),!-O,!=O)'),
    ('cyano','C(#N)'),
    ('n_cyano','N(#C)'),
    ('nitro','C(-N(-O(!-C,!-N,!-H),=O))'),
    ('amido','C(=O,-N(!=C))'),
    ('imine','C(=N(=C,!-O),!-O,!=O)'),
    ('methyl','C(-H,-H,-H)'),
    ('ethyl','C(-H,-H,-methyl)'),
    ('propyl','C(-H,-H,-C(-H,-H,-methyl))'),
    ('isopropyl','C(-H,-methyl,-methyl)'),
    ('butyl','C(-H,-H,-C(-H,-H,-C(-H,-H,-methyl)))'),
    ('isobutyl','C(-H,-C(-H,-H,-methyl),-methyl)'),
    ('tert_butyl','C(-methyl,-methyl,-methyl)'),
    ('methylene','C(-H,-H,!-H)'),
    ('thioether','C(-S(-C(!=O)),!=O)'),
    ('thioester','C(=O,-S(-C))'),
    ('reverse_thioester','C(-S(-C(=O)))'),
    ('primary_thiol','C(-S(-H),-H,-H)'),
    ('secondary_thiol','C(-S(-H),-H,!-H)'),
    ('tertiary_thiol','C(-S(-H),!-H,!-H)'),
    ('sulfonamide','C(-N(-S(=O,=O,-C)))'),
    ('thiobridge','C(-S(-S(-C)))'),
    ('methylthioether','C(-S(-C,-methyl))'),
    ('phosphoester','C(-O(-P(=O,-O,-O)),!=O,!-N)'),
    ('hydroxy_benzene','C9(=C(-C(=C(-C(=C(-C9))))),-O(-H))'),
    ('phenolic_ester','C9(=C(-C(=C(-C(=C(-C9))))),-O(-C(=O)))'),
    ('phenolic_phosphate_ester','C9(=C(-C(=C(-C(=C(-C9))))),-O(-P(=O)))'),
    ('o-glycoside','C(-O(-C(-O(-C(!~O)),!-O)))'),
    ('n-glycoside','C(-O(-C(-N(-C(!~O)),!-O)))'),
    ('benzene','C9(-C(=C(-C(=C(-C(=C9))))))')    
    )

FUNCTIONAL_GROUP_TOPOLOGIES = [Pattern(n, p) for n, p in FUNCTIONAL_GROUP_PATTERNS]
ALPHA_TOPOLOGIES = [Pattern(n, p) for n, p in ALPHA_GROUP_PATTERNS]
