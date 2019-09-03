#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    for l in open(args.file):
        if l.startswith('# ') or l.strip() == '//':
            print(l)
        if l.startswith('#=GF'):
            print(l)
        elif '|' in l and '#' not in l:
            #0000|Chr4_88684848__88684954             GUGCUCg........
            # print(l)
            # print(l.split('|'))
            x = l.split('|')[1]
            print(x.strip())
