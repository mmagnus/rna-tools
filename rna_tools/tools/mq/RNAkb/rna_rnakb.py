#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing NAST potential via Pyro.
"""
from __future__ import print_function

import argparse
import os
from rna_tools.tools.mq.RNAkb.RNAkb import RNAkb


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
        wrapper = RNAkb(sandbox=True)
        result = wrapper.run(f, '5pt', args.verbose)
        print(f + ',' + str(result))
        wrapper.cleanup()
