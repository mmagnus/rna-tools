#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_pairs_diff.py - get a diff of pairs

Usage::

    $ rna_pairs_diff.py pistol_dca_all.bp pistol.bp
    # of ec_paris: 31
    # of ssbps   : 18
    dalta#       : 13
    [[4, 32], [6, 9], [6, 36], [6, 39], [9, 39], [13, 32], [16, 17], [17, 18], [22, 49], [29, 58]]

"""
from __future__ import print_function

import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('pairs1', help="a list of pairs, A")
    parser.add_argument('pairs2', help="a list of pairs to subtract, A-B, results in C"
                        "(all pairs that are in A and are not in B")
    parser.add_argument("-v", "--verbose", action='count',
                        help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    pairs1 = eval(open(args.pairs1).read().strip())
    pairs2 = eval(open(args.pairs2).read().strip())

    # go ever dca and keep if not in ss
    pairsnonss = []
    for pair in pairs1:
        if pair in pairs2:
            pass
        else:
            pairsnonss.append(pair)
    print('# of ec_paris:', len(pairs1))
    print('# of ssbps   :', len(pairs2))
    print('dalta#       :', len(pairsnonss))
    pairsnonss.sort()
    print(pairsnonss)
