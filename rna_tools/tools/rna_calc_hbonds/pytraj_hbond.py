#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        print(f)
        import pytraj as pt
        pdb  = pt.load(f)
#          print(pdb)
        hbonds = pt.search_hbonds(pdb)
        ic(hbonds)
        for h in hbonds:
            print(h)
        print(hbonds)
