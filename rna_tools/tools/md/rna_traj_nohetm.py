#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
hack::

     parallel rna_traj_nohetm.py {} ::: *MD.*

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
        npdb = ''
        for l in open(f):
            if l.startswith('HETATM'):
                pass
            else:
                npdb += l
            if l.startswith('MODEL'):
                print(l.strip())

        nf = f.replace('.pdb', '_rmHtm.pdb')
        with open(nf, 'w') as f:
            f.write(npdb)
        print('save %s' % nf)
