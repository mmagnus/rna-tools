#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_restrs2mqapRNA.py - convert pairs to mqapRNA restraints

Example::

    $ rna_restrs2mqapRNA.py dca_pairs.txt
    d:A1 - A86 < 7
    d:A2 - A85 < 7
    d:A3 - A83 < 7
    d:A3 - A84 < 7
    d:A4 - A49 < 7

"""
from __future__ import print_function
import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('pairs', help="a file with [[2, 172], [3, 169], [12, 32], [13, 31]]")
    parser.add_argument("--offset", help="can be -10", default=0, type=int)
    parser.add_argument("--weight", default=3, type=float, help="weight")
    parser.add_argument("--dist", default=7, help="distances, for MOHCA use 25", type=int)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    pairs = eval(open(args.pairs).read().strip())
    if args.verbose:
        print('# of pairs:', len(pairs))

    for pair in pairs:

        a, b = pair
        a += args.offset
        b += args.offset

        a = str(a)
        b = str(b)

        print('d:A' + a + '-A' + b + ' < 7 ')
