#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_pairs2SimRNArestrs.py - convert pairs to SimRNA restraints

Example::

    $ rna_paris2rfilter.py rp06_MohPairs.bp
    d A5-A42< 100 1
    d A5-A42< 100 1
    d A11-A26< 100 1
    d A11-A26< 100 1

Input::

    $ head rp06_MohPairs.bp
    [[5, 42], [11, 26], [24, 89], [26, 44], [26, 75],

"""
from __future__ import print_function
import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('pairs', help="a file with [[2, 172], [3, 169], [12, 32], [13, 31]]")
    parser.add_argument("--offset", help="can be -10", default=0, type=int)
    parser.add_argument("--weight", default=3, type=float, help="weight")
    parser.add_argument("--distance", default=10, type=float, help="distance")
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

        # hmmm.. why 2 lines here?
        print('d:A' + a + '-A' + b + ' < %.2f %.2f' % (args.distance, args.weight))
        #print('d A' + a + '-A' + b + ' < %f %f' % (args.distance, args.weight)
