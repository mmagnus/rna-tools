#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lazy shortcut for::

    rna_extract.static.macosclangrelease -in::file::silent *out # be default

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
    #parser.add_argument('-', "--", help="", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--file", help="", default="*.out")#, nargs='+')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    import os
    os.system('rna_extract.static.macosclangrelease -in::file::silent ' + args.file) # *.out
