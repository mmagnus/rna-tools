#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_ss_get_bps.py - get a list of base pairs for a given "fasta ss" file.

Input file::

    cat ade_ss.fa
    >1y26
    CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
    (((((((((...((((((.........))))))........((((((.......))))))..)))))))))%

Usage::

    $ rna_ss_get_bps.py ade_ss.fa --offset 12
    [[13, 83], [14, 82], [15, 81], [16, 80], [17, 79], [18, 78], [19, 77], [20, 76], [21, 75], [25, 45], [26, 44], [27, 43], [28, 42], [29, 41], [30, 40], [54, 72], [55, 71], [56, 70], [57, 69], [58, 68], [59, 67]]

Now it also work with pseudoknots.
"""
from __future__ import print_function

import argparse
from rna_tools.SecondaryStructure import parse_vienna_to_pairs


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file', help="file in the Fasta format")
    parser.add_argument('--offset', help="offset", type=int)
    parser.add_argument("-v", "--verbose", action='count',
                        help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    id, seq, ss = open(args.file).read().strip().split('\n')
    pairs, pairs_pk = parse_vienna_to_pairs(ss)

    if args.offset:
        npairs = []
        for pair in pairs:
            npairs.append([pair[0] + args.offset, pair[1] + args.offset])
        for pair in pairs_pk:
            npairs.append([pair[0] + args.offset, pair[1] + args.offset])
    else:
        npairs = pairs + pairs_pk
    npairs.sort()
    print(npairs)
