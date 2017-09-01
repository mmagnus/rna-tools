#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_ex2x.py - analyze an evolutionary coupling file.

Files can be downloaded from https://marks.hms.harvard.edu/ev_rna/, e.g. RF00167.EC.interaction.csv

``--pairs``::

    $ rna_ex2x.py RF00167.EC.interaction_LbyN.csv --pairs
    [18, 78],[31, 39],[21, 75],[30, 40],[28, 42],[27, 43],[59, 67],[54, 72],[57, 69],[25, 45],[29, 41],[17, 79],[26, 44],[16, 80],[14, 82],[19, 77],[55, 71],[
        15, 81],[34, 63],[56, 70],[58, 68],[35, 63],[26, 45],[35, 64],[32, 39],[54, 73],[24, 74],[16, 82],[24, 45],[24, 43],[32, 36],[25, 48],[48, 82],[36, 48],
"""
from __future__ import print_function

import csv
import pandas
import sys
import argparse


def get_parser():
    """Get parser."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--sep", default=r",", help="separator")
    parser.add_argument("--chain", default="A", help="chain")
    # parser.add_argument("--offset", default=0)
    parser.add_argument('interaction_fn', help='interaction file')
    parser.add_argument('--ec-pairs',  action='store_true')
    parser.add_argument('--ss-pairs',  help="file with secondary structure base pairs")
    parser.add_argument('--pairs-delta', help="delta: ec-bp - ss-paris", action="store_true")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    df = pandas.read_csv(args.interaction_fn, sep=",")
    pairs = []
    for index, row in df.iterrows():
        if str(row['pdb_resid1']) != 'nan':
            # ? uggly hack to 78.0 -> 78
            if not args.ec_pairs:
                print('d:' + args.chain + str(row['pdb_resid1']).replace('.0', '') +
                      '-' + args.chain + str(row['pdb_resid2']).replace('.0', '') + ' < 100 1')
            else:
                pairs.append([int(row['pdb_resid1']), int(row['pdb_resid2'])])
                pairs.sort()

                # this also worked, but it was an ugly print
                # print('[' + str(row['pdb_resid1']).replace('.0', '') + ', ' +
                #    str(row['pdb_resid2']).replace('.0', '') + ']', end=',')

    if args.ec_pairs:
        print('pairs:')
        print(pairs)

    if args.ss_pairs:
        ssbps = eval(open(args.ss_pairs).read().strip())
        ssbps.sort()
        print('ssbps:')
        print(ssbps)

    if args.pairs_delta:
        # go ever dca and keep if not in ss
        pairsnonss = []
        for pair in pairs:
            if pair in ssbps:
                pass
            else:
                pairsnonss.append(pair)
        print('# of ec_paris:', len(pairs))
        print('# of ssbps   :', len(ssbps))
        print('dalta#       :', len(pairsnonss))
        print(pairsnonss)
