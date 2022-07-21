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
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-o", "--output", help="", default="_merged_.csv")
    parser.add_argument('--sep', help="default is ,; can be also '\t'", default=",")
    parser.add_argument('mergeon', help="merge on column", default="")
    parser.add_argument("files", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    files = args.files
    f1 = files.pop()
    merged = pd.read_csv(f1, delimiter=args.sep)
    for f in files:
        dft = pd.read_csv(f, delimiter=args.sep)
        merged = pd.merge(merged, dft, on=args.mergeon)#, how=args.how, validate=args.validate)
    print(merged)
    print('saved ', args.output)
    merged.to_csv(args.output, index=False)

