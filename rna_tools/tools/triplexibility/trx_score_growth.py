#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    # parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    dat = """aua,0
    auc,0
    auG,0
    aua,0.8
    Cca,0
    CGa,0
    CcG,0
    Ccu,0
    CGa,0
    CGc,0
    CGG,1
    CGu,1
    Cua,0
    Cuc,0
    Cug,0
    Cuu,0.5"""

    for l in dat.split('\n'):
        s, g = l.split(',')
        s = s.strip()
        # load data
        d = open('/Users/magnus/Desktop/trx-farfar/t1-2/' + s + '/' + s + '_short.csv').read() # _full
        print(s + ',' + g + ',' + d.replace('source,target,n\n', '').replace('\n', ',')[:-1]) # to remove last coma
