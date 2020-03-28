#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
import subprocess
import os
from rna_tools.Seq import RNASequence


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--method', help="rnaeval", default="rnaeval")#, required=True)
    parser.add_argument("--seq", help="", default="")
    parser.add_argument("--file", help="", default="", nargs='+')
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

def hr(t):
    print
    if t:
        print(t)
    print('-' * 80)
    
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    cst = args.cst

    for file in args.file:

        hr(file)

        if args.seq:
            seq = args.seq

        if args.file:
            f = open(file)
            header = f.readline()
            seq_string = f.readline().strip()
            if args.cstinfile:
                cst = f.readline().strip()
                assert cst, 'Cst line empty, check the input file'

        try:
            print(seq_string)
        except NameError:
            parser.print_usage()
            parser.exit(1)

        if args.cst:
            cst = args.cst
        if cst:
            print(cst, '<= cst')

        if args.method == 'rnaeval' or args.method == 'all':
            seq = RNASequence(seq_string)
            energy = seq.eval(cst, verbose=args.verbose)
            print('rnaeval:', energy)

        if args.method == 'mcfold' or args.method == 'all':
            seq = RNASequence(seq_string)
            if cst:
                en, ss, comment = seq.predict_ss('mcfold', constraints=cst, verbose=args.verbose)
            else:
                en, ss, comment = seq.predict_ss('mcfold', verbose=args.verbose)
            #if not ss:
            #    ss = 'x'
            print('mcfold:', en, comment)
