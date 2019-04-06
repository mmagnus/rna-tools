#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import pandas as pd
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('csv1', help="", default="")
    parser.add_argument('--sep1', help="default is ,; can be also '\t'", default=",")
    parser.add_argument('csv2', help="", default="")
    parser.add_argument('--sep2', help="default is ,; can be also '\t'", default=",")
    parser.add_argument('mergeon', help="merge on column", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('--output', help="output csv", default="merged.csv")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    df1 = pd.read_csv(args.csv1, delimiter=args.sep1)
    df2 = pd.read_csv(args.csv2, delimiter=args.sep2)
    print(df1)
    print(df2)

    merged = pd.merge(df1, df2, on=args.mergeon)
    print(merged)
    merged.to_csv(args.output, index=False)
    print('Saved to %s' % args.output)
