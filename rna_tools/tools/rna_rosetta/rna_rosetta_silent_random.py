#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
import os
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", '--output', help="output file", default="")
    parser.add_argument("--keep-order",
                        action="store_true", help="keep order of S* like in the input file")
    parser.add_argument("n", help="", default=10)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    n = int(args.n)

    ic.disable()
    if args.verbose:
        ic.enable()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        fn = f
        print(f)
        f = open(f)
        seq = f.readline()
        headers = f.readline()
        remark = f.readline()
        strucs = []
        struc = ''
        ids = []

        for line in f:
            if len(line) > 5 and line[:5] == 'SCORE':  # start new S
                if struc:
                    strucs.append(struc)
                struc = line # set line
                id = line.split()[-1]
                ids.append(id)
                ic(id)
            else:
                struc += line
        strucs.append(struc) # get the last one

        # keep order
        from numpy import random
        selected = random.choice(range(0, len(ids)), size=n, replace=False)
        ic(selected, ids)
        sall = ''
        hall = ''
        if args.keep_order:
            for i, s in enumerate(strucs):
                if i in selected:
                    ic(i, selected)
                    sall += strucs[i]
                    hall += ids[i] + ' '
        else:
            for sel in selected:
                sall += strucs[sel]
                hall += ids[sel] + ' '

        if args.output:
            nfn = args.output
        else:
            nfn = fn + '.rand' + str(n) + '.out'
        with open(nfn, 'w') as f:
            print('Random of %i into %s: %s' % (n, nfn, hall))
            f.write(seq)
            f.write(headers)
            f.write(remark)
            f.write(sall)
