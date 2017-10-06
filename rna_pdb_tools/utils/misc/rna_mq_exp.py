#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import pandas as pd
import argparse
import os


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('csv', help="", default="")
    parser.add_argument('--sep', help="default is ,; can be also '\t'", default=",")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('--output', help="output csv", default="merged.csv")
    parser.add_argument('--inf', help="inf", default="inf.csv")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    header = 'fn,mqapRNA,rst,'
    if args.inf:
        header += 'inf,mqapRNArst'
    else:
        header += 'mqapRNArst'

    print(header)

    for l in open(args.csv):
        if l.startswith('fn'):
            pass
        else:
            # otherwise there is a problem with float(inf)
            inf = 0

            fn = l.split(',')[0].strip()
            # parser inf.csv file to get inf score
            # search for a file in inf.csv with the same fn as in csv
            if os.path.isfile(args.inf):
                for r in open(args.inf):
                    # target.pdb.outCR,2gdi_tpp-022e50b3-thrs8.00A_clust04X.pdb.outCR,0.848,0.000,0.941,0.000,1.000,0.885,0.000,0.000
                    cells = r.split(',')
                    if not r.startswith('target,fn'):
                        inf_fn = cells[1].replace('.outCR', '')
                        if fn == inf_fn :
                            inf = 1 - float(cells[4]) # inf_WC # 1 -> expressed as a penalty!
                            break

            mq = float(l.split(',')[1].strip())
            exp = 1 - float(l.split(',')[2].strip()) # as penalty right now

            if args.inf:
                mq_exp = float(mq) + exp + float(inf)
            else:
                mq_exp = float(mq) + exp #
            # +(1 * ((1 - exp) * mq))
            # #mq_exp = float(mq) +(1 * ((1 - exp) * mq))
            if args.inf:
                print(','.join([fn, str(mq), str(exp), str(inf), str(mq_exp)]))
            else:
                print(','.join([fn, str(mq), str(exp), str(mq_exp)]))
