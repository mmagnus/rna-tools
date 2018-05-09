#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
from moderna import *
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    ## parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    # mutate
    # --mutate 'A:1C+A:10A' # ok, lets keep plus!
    # load each chain and mutate etc. and merge them into one
    ## m = load_model('input/1xjr_A1-4.pdb')
    ## exchange_single_base(m['1'],'C')
    ## exchange_single_base(m['2'],'C')
    ## m.write_pdb_file('output/1xjr_A1-4_A1C.pdb')

    ## m = load_model('input/CG.pdb')
    ## exchange_single_base(m['9'],'G')
    ## exchange_single_base(m['41'],'C')
    ## m.write_pdb_file('output/GC2CG.pdb')

    # this is ok
    ## m = load_model('input/205d_rmH2o_1chain.pdb')
    ## exchange_single_base(m['4'],'G')
    ## exchange_single_base(m['21'],'C')
    ## m.write_pdb_file('output/205d_rmH2o_1chain_mutant.pdb')


    m = load_model('input/205d_rmH2o.pdb')
    exchange_single_base(m['4'],'G')
    m.write_pdb_file('output/205d_rmH2o_mutant_A.pdb')

    m = load_model('input/205d_rmH2o.pdb', 'B')
    exchange_single_base(m['21'],'C')
    m.write_pdb_file('output/205d_rmH2o_mutant_B.pdb')

    # _A.pdb
    # _B.pdb
    # and then
    # test what if you GC -> CG
