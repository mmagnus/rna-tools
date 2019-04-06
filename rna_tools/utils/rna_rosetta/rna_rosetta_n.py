#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_roseta_n.py - show me # of structure in a silent file

Example::

     $ rna_rosetta_n.py ade.out
     21594

"""
from __future__ import print_function
import subprocess
import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('file', help='ade.out')
    return parser


def get_no_structures(file):
    p = subprocess.Popen('cat ' + file + ' | grep SCORE | wc -l', shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stderr = p.stderr.read().strip()
    if stderr:
        print(stderr)
    return int(p.stdout.read().strip()) - 1


def run():
    """Pipline for modeling RNA."""
    args = get_parser().parse_args()
    ns = get_no_structures(args.file)
    if args.verbose:
        print(args.file.ljust(30), end='')
    print(ns)


# main
if __name__ == '__main__':
    run()
