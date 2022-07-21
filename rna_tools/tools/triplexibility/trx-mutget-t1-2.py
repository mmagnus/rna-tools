#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    aa='x,c,c,a,c'
    bb='x,g,u,u,c'
    cc='g,u,a,c'

    csv = open(args.file).read().split('\n')

    for a, b in zip(aa.split(','), bb.split(',')):
         if a != 'x':
             for c in cc.split(','):
                 seq = 'A:22' + a + ',C:60' + b + ',B:52' + c # search for the last column
                 cmd = "rna_pdb_toolsx.py --mutate '" + seq + "' " + args.file + " > " + "" + a + b + c + '-' + args.file + seq.replace(':', '-').replace(',', '-') + "-rpr.pdb"
                 print(seq)
                 print(cmd)
                 import os
                 os.system(cmd)
