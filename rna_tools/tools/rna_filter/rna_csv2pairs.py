#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_csv2pairs.py - get pairs based on a csv file.

Usage::

    $ rna_csv2pairs.py rp06_MohPairs.txt
    [[5, 42], [11, 26], [24, 89], [26, 44], [26, 75], [35, 75],

Input::

    $ head rp06_MohPairs.txt
    11,26
    26,44
    5,42

"""
from __future__ import print_function

import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('fn', help="csv file", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    pairs = []
    for l in open(args.fn):
        if not l.startswith('#'):
            a, b = l.strip().split(',')
            a = int(a)
            b = int(b)
            pairs.append([a, b])
    pairs.sort()
    print(pairs)
