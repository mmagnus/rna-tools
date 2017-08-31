#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
`--pairs``::

    $ rna_csv2restrs.py RF00167.EC.interaction_LbyN.csv --pairs
    [18, 78],[31, 39],[21, 75],[30, 40],[28, 42],[27, 43],[59, 67],[54, 72],[57, 69],[25, 45],[29, 41],[17, 79],[26, 44],[16, 80],[14, 82],[19, 77],[55, 71],[
        15, 81],[34, 63],[56, 70],[58, 68],[35, 63],[26, 45],[35, 64],[32, 39],[54, 73],[24, 74],[16, 82],[24, 45],[24, 43],[32, 36],[25, 48],[48, 82],[36, 48],

old print::

    rna_EC2restrs.py RF00167.EC.interaction_LbyN.csv --pairs
    [18, 78],[31, 39],[21, 75],[30, 40],[28, 42],[27, 43],[59, 67],[54, 72],[57, 69],[25, 45],[29, 41],[17, 79],[26, 44],[16, 80],[14, 82],[19, 77],[55, 71],[
        15, 81],[34, 63],[56, 70],[58, 68],[35, 63],[26, 45],[35, 64],[32, 39],[54, 73],[24, 74],[16, 82],[24, 45],[24, 43],[32, 36],[25, 48],[48, 82],[36, 48],

vs new list (thanks to list, it's sorted now!)::

    rna_EC2restrs.py RF00167.EC.interaction_LbyN.csv --pairs
    [[14, 82], [15, 81], [16, 80], [16, 82], [17, 79], [18, 78], [19, 77], [21, 75], [24, 43], [24, 45], [24, 74], [25, 45], [25, 48], [26, 44], [26, 45], [27, 43], [28, 42], [
        29, 41], [30, 40], [31, 39], [32, 36], [32, 39], [34, 63], [35, 63], [35, 64], [36, 48], [48, 82], [54, 72], [54, 73], [55, 71], [56, 70], [57, 69], [58, 68], [59, 67]]

"""
from __future__ import print_function

import csv
import pandas
import sys
import argparse


def get_parser():
    parser = argparse.ArgumentParser()  # usage="prog [<options>] <pdb files: test_data/*>")
    parser.add_argument("--sep", default=r",")
    parser.add_argument("--chain", default="A")
    # parser.add_argument("--offset", default=0)
    parser.add_argument('interaction_fn')
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
