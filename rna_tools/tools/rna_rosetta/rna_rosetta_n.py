#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_roseta_n.py - show me # of structure in a silent file

Example::

    $ rna_rosetta_n.py ade.out
    21594

    $ rna_rosetta_n.py *out
    default.out 100
    default1.out 200
    default2.out 200
    default3.out 100

"""
from __future__ import print_function
import subprocess
import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-0', '--zero', action='store_true')
    parser.add_argument('file', help='ade.out', nargs='+')
    return parser


def get_no_structures(file):
    p = subprocess.Popen('cat ' + file + ' | grep SCORE | wc -l', shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stderr = p.stderr.read().strip()
    if stderr:
        print(stderr)
    return int(p.stdout.read().strip()) - 1


def run(f):
    """Pipline for modeling RNA."""
    args = get_parser().parse_args()
    import os
    if os.path.isfile(f):
        ns = get_no_structures(f)
        if args.verbose:
        #    print(args.file.ljust(30), end='')
             print(f, ns)
        else:
             print(ns)        
    else:
        if args.zero:
            print(0)
        else:
            print(f, 'missing')
        
# main
if __name__ == '__main__':


    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        run(f)
