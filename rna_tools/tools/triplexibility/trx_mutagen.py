#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("ref", help="", default="")
    parser.add_argument("--fig", help="", default="fig.py")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    exec(open(args.fig).read())

    for a, b in zip(aa.split(','), bb.split(',')):
         if a != 'x':
             for c in cc.split(','):
                 seqm = a + b + c
                 #seq = 'A:1' + a + '+2' + b + '+3' + c # search for the last column
                 seq = 'A:1' + a + '+2' + b + '+3' + c # search for the last column
                 f2 = args.ref + "_" + a + b + c + ".pdb"
                 cmd = "rna_pdb_tools.py --mutate '" + seq + "' " + args.ref + " > " + f2
                 print(seq)
                 print(cmd)
                 os.system(cmd)
                 continue
             
                 f3 = f2.replace('.pdb', '_rpr.pdb')
                 cmd = "rna_pdb_tools.py --rpr " + f2 + " > " + f3
                 print(cmd)
                 os.system(cmd)
                 # --files ../../db/triples-all-v2-rpr/*
                 csvf = a + b + c + ".csv"
                 cmd = "triplexibility2.py -t " + f3 + " --save --sort --result " + csvf
                 print(cmd)
                 os.system(cmd)
