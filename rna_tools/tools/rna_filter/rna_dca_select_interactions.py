#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""select top interactions

I implemented correction to sync these results with excel files from the Marks lab.
"""
from __future__ import print_function
import argparse
import pandas as pd
import sys


def sortdata(filename):
    df = pd.read_csv(filename, sep=" ", names=['i', ' ', 'j', ' ', ' ', 'scores'])
    df2 = df.sort_values(by=['scores'], ascending=False)
    maxvalue = max(df2[['j']].values)
    l = int(maxvalue)
    print('l', l)
    L = l / 2
    #L = l
    print('L', L)
    df3 = df2[0:L]
    df4 = df3[['i', 'j', 'scores']]
    # correction for -1 ?
    df4['i'] -= 1
    df4['j'] -= 1

    print(df4)
    csvfn = filename + "_" + str(L) + "_selected.csv"
    print('Output file created: ' + csvfn)
    df4.to_csv(csvfn, sep=" ", index=False)


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file', help="plm scores, direct output from plm", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    sortdata(args.file)
