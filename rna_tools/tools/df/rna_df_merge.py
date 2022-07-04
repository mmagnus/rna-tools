#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
import pandas as pd
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
    parser.add_argument("--drop-col", help="", default="")
    parser.add_argument("-o", "--output", help="", default="_merged_.csv")
    parser.add_argument("files", help="", default="", nargs='+')
    
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    args.sep = ','
    dfs = []
    for f in args.files:
        df = pd.read_csv(f, delimiter=args.sep)
        if args.drop_col:
            df.drop(args.drop_col, axis=1, inplace=True)
        print(df)
        dfs.append(df)
    dfs = pd.concat(dfs, axis=1)
    print(dfs)
    print('saved ', args.output)
    dfs.to_csv(args.output, index=False)
