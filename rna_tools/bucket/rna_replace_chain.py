#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
from rna_tools_lib import RNAStructure
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser




if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    s = RNAStructure('input/205d_rmH2o.pdb')
    chain_id = 'A'
    insert = RNAStructure('output/205d_rmH2o_mutant_A.pdb')

    out = replace_chain('input/205d_rmH2o.pdb',
                        'output/205d_rmH2o_mutant_A.pdb',
                            chain_id)
    print(out)
