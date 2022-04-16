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
    parser.add_argument("-o", '--output-dir', help="", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
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
        strucs.append(struc)

        for i, s in zip(ids, strucs):
            if args.output_dir:
                fn = args.output_dir + os.sep + str(i) + '.out'
            else:
                fn = str(i) + '.out'
            with open(fn, 'w') as f:
                print('Split into %s' % fn)
                f.write(seq)
                f.write(headers)
                f.write(remark)
                f.write(s)

