#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import os
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    txt = """cSW_cSH_CAA
cSW_cSH_CAG
cSW_tHW_CAA
cSW_tHW_CAC
cSW_tSW_UGC
cSW_tSS_CAA
cSW_tSS_CAG
cWH_cSH_UGC
cWW_cSH_UGC
cWW_tHS_CAA
cWW_tHS_CAC
cWW_tSW_CAC
tHH_cSH_CAC
tHH_tSS_CAC
tHH_tSS_GCA
tSH_tSW_CAC
tSW_cSH_UGC
tSW_tHS_CAA
tSW_tHS_CAC
tSW_tHS_CUG
tWH_cSH_CAA
tWH_cSH_CAG
tWH_tSW_UGC
tWW_cSH_CAC
tWW_tHW_CAA
tWW_tHW_CAC
tWW_tHW_CUG
tWW_tSW_UGC
tWW_tSS_CAA
tWW_tSS_CAG
"""
    txt = """
tHH_tSS_GCA
"""
    for t in txt.split('\n'):
        print(t)
        os.system('wget http://rna.bgsu.edu/triples/Data/v1.4/Models/Triple_' + t + '.pdb')


