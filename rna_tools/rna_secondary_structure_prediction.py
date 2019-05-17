#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
from Seq import RNASequence
import argparse
import sys
import subprocess
import os
reload(sys)
sys.setdefaultencoding('utf8')

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--method', help="cyclefold, mcfold", default="cyclefold")#, required=True)
    parser.add_argument("--seq", help="", default="")
    parser.add_argument("--file", help="", default="")
    parser.add_argument("--cstinfile", action="store_true", help="the 3rd file should have cst line, e.g. '((....))'")
    parser.add_argument("--cst", help="--cst '((....))'")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


def exe2(cmd):
    os.system(cmd)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    cst = args.cst

    if args.seq:
        seq = args.seq

    if args.file:
        f = open(args.file)
        header = f.readline()
        seq = f.readline().strip()
        if args.cstinfile:
            cst = f.readline().strip()
            assert cst, 'Cst line empty, check the input file'

    try:
        print(seq)
    except NameError:
        parser.print_usage()
        parser.exit(1)

    if args.cst:
        cst = args.cst
    if cst:
        print(cst, '<= cst')

    if args.method == 'cyclefold':
        seq = RNASequence(seq)
        ss, en = seq.predict_ss('rnastructure_CycleFold', verbose=args.verbose)
        print(en, ss)

    if args.method == 'mcfold':
        seq = RNASequence(seq)
        if cst:
            en, ss = seq.predict_ss('mcfold', constraints=cst, verbose=args.verbose)
        else:
            en, ss = seq.predict_ss('mcfold', verbose=args.verbose)
        #if not ss:
        #    ss = 'x'
        print(ss, en)
