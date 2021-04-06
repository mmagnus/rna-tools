#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
import glob

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-d", "--debug",
                        action="store_true", help="be verbose")
    parser.add_argument("seq", help="")#, nargs='+')
    parser.add_argument("--edge", help="e.g., cWW_cHS")
    return parser


def hr(t):
    l = len(t)
    half = int(l/2)
    print('-' * (40 - half) + ' ' + t + ' ' + '-' * (40 - half))

def fopen(filename, verbose=0):
    import subprocess
    cmd = 'mdfind -name "' + filename + '"'
    if verbose:
        print(cmd)
    out = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).stdout.read().decode()
    first_hit = out.split('\n')[0]
    return open(first_hit)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    edge = args.edge
    seq = args.seq
    
    import pandas as pd
    df = pd.read_csv("../triples/triple-db.csv")
    triple = edge + '_' + seq

    x = df[df['triple'] == triple]
    r3 = fopen('Triple_' + edge + '_' + seq  + '_rpr.pdb.3dcnn.csv').read()
    print(triple + ':', round(int(x['instances']) / 18, 2), r3.split(',')[1].strip())


    



