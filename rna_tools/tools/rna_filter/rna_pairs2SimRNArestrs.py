#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_pairs2SimRNArestrs.py - convert pairs to SimRNA restraints

Example::

    $ rna_pairs2SimRNArestrs.py rp06_pairs_delta.txt -v
    # of pairs: 42
    SLOPE A/2/MB A/172/MB 0 6 1
    SLOPE A/2/MB A/172/MB 0 7 -1
    SLOPE A/3/MB A/169/MB 0 6 1
    SLOPE A/3/MB A/169/MB 0 7 -1
    SLOPE A/12/MB A/32/MB 0 6 1

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
    parser.add_argument("--well", help="well instead of slope", action="store_true", default=False)
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

        if args.well:
            print('WELL A/' + a + '/MB A/' + b + '/MB 0 ' + str(args.dist) + ' ' + str(args.weight))
        else:
            print('SLOPE A/' + a + '/MB A/' + b + '/MB 0 ' + str(args.dist - 1) + ' ' + str(args.weight))
            print('SLOPE A/' + a + '/MB A/' + b + '/MB 0 ' + str(args.dist) + ' -' + str(args.weight))
