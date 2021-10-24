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
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        import mdtraj as md
        print(f, '---------------------')
        t = md.load(f)
        print(md.kabsch_sander(t))
        hb = md.wernet_nilsson(t)
        print(hb)
        print(len(hb[0]))        
        hb = md.baker_hubbard(t)
        print(hb)
        print(len(hb))
