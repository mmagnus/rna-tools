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
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    print('fn,mqapRNA,rst,mqapRNArst')
    for l in open(args.csv):
        if l.startswith('fn'):
            pass
        else:
            fn = l.split(',')[0].strip()
            mq = float(l.split(',')[1].strip())
            exp = float(l.split(',')[2].strip())
            mq_exp = float(mq) +(1 * ((1 - exp) * mq))
            print(','.join([fn, str(mq), str(exp), str(mq_exp)]))
