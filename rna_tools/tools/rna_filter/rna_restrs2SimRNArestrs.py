#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Keep only for reference. TO BE REMOVE!!!!

SLOPE A/23/MB C/45/MB 0 6 2.0
SLOPE A/23/MB C/45/MB 0 7 -2.0

A16-A82 < 100 1

"""
from __future__ import print_function
import sys


import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--offset", help="can be -10", default=0, type=int)
    parser.add_argument("--weight", help="", default=1, type=float)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('file')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    # 13 -> 1
    offset = args.offset
    fn = args.file
    weight = args.weight
    #
    for l in open(fn):
        if l.startswith('#'):
            continue

        # d:A58-A68 < 100 1
        l = l.replace('d:', '')
        l = l.replace('A', '')
        pair, rest = l.split('<')
        a, b = pair.split('-')

        a = int(a.strip())
        a += offset
        b = int(b.strip())
        b += offset

        a = str(a)
        b = str(b)

        print('SLOPE A/' + a + '/MB A/' + b + '/MB 0 6 ' + str(weight))
        print('SLOPE A/' + a + '/MB A/' + b + '/MB 0 7 -' + str(weight))
