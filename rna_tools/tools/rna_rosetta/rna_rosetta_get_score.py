#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

SEQUENCE: gaugaucagcaguuccccugcauaaggaugaaccggcuuagauc
SCORE:     score     fa_atr     fa_rep    fa_intra_rep    lk_nonpolar    fa_elec_rna_phos_phos    rna_torsion    suiteness_bonus    rna_sugar_close    fa_stack    stack_elec    geom_sol_fast    hbond_sr_bb_sc    hbond_lr_bb_sc    hbond_sc        ref    free_suite    free_2HOprime    intermol    other_pose    missing    loop_close    linear_chainbreak description
REMARK BINARY SILENTFILE    FULL_MODEL_PARAMETERS  FULL_SEQUENCE gaugaucagcaguuccccugcauaaggaugaaccggcuuagauc  CONVENTIONAL_RES_CHAIN A:52-85,B:20-29  CUTPOINT_OPEN 34  DOCK_DOMAIN 1-44  FIXED_DOMAIN 1-5,7,9-37,39-44  INPUT_DOMAIN 1-5,7,9-37,39-44  SAMPLE 6,8,38  WORKING 1-44
SCORE:   171.254   -128.326     58.303           1.812          2.997                   18.602        118.895              0.000             36.107     -81.386        -4.746           20.132             0.000            -3.060     -30.102    153.360        -2.000            0.000       2.300         0.000          1         7.878                0.488    S_000001
RES_NUM A:52-56 A:58-85 B:20-29 S_000001
FOLD_TREE  EDGE 1 4 -1  EDGE 4 5 -1  JEDGE 1 28 2 C4 C2  INTRA_RES_STUB  JEDGE 4 42 3 C4 C2  INTRA_RES_STUB  EDGE 28 10 -1  EDGE 28 33 -1  EDGE 42 38 -1  EDGE 42 43 -1  JEDGE 10 35 1 C4 C4  INTRA_RES_STUB  JEDGE 43 6 4 C2 C2

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    c = 2
    scores = []
    for f in args.file:
        row = 'id,score,fn'
        print(row)
        row = '1,'            
        for l in open(f):
            if l.startswith('SCORE:     score'):
                continue
            if l.startswith('SCORE'):
                score = float(l.split()[1])
                #ic(score)
                #print(score)
                row += str(score)
                scores.append(score)

            if 1: # False:
                if l.startswith('RES_NUM'):
                    fn = l.split()[-1] + '.pdb'
                    #ic(fn)
                    row += ',' + fn
                    print(row)
                    row = str(c) + ','
                    c += 1
    scores.sort()
    print(scores)
            
