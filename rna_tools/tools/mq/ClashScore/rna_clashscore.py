#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing NAST potential via Pyro.
"""
from __future__ import print_function

import argparse
import os
from rna_tools.tools.mq.ClashScore.ClashScore import ClashScore


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
        wrapper = ClashScore()
        result = wrapper.run(f, args.verbose)
        print(f + ', ' + str(result))
        wrapper.cleanup()
