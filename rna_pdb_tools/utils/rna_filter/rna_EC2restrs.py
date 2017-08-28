#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
`--pairs``::

    $ rna_csv2restrs.py RF00167.EC.interaction_LbyN.csv --pairs
    [18, 78],[31, 39],[21, 75],[30, 40],[28, 42],[27, 43],[59, 67],[54, 72],[57, 69],[25, 45],[29, 41],[17, 79],[26, 44],[16, 80],[14, 82],[19, 77],[55, 71],[15, 81],[34, 63],[56, 70],[58, 68],[35, 63],[26, 45],[35, 64],[32, 39],[54, 73],[24, 74],[16, 82],[24, 45],[24, 43],[32, 36],[25, 48],[48, 82],[36, 48],

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
    parser.add_argument('--pairs',  action='store_true')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    df = pandas.read_csv(args.interaction_fn, sep=",")
    for index, row in df.iterrows():
        if str(row['pdb_resid1']) != 'nan':
            # ? uggly hack to 78.0 -> 78
            if not args.pairs:
                print('d:' + args.chain + str(row['pdb_resid1']).replace('.0', '') +
                      '-' + args.chain + str(row['pdb_resid2']).replace('.0', '') + ' < 100 1')
            else:
                print('[' + str(row['pdb_resid1']).replace('.0', '') + ', ' +
                      str(row['pdb_resid2']).replace('.0', '') + ']', end=',')
