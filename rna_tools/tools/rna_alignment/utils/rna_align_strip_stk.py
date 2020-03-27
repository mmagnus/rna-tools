#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    for l in open(args.file):
        # AF035635.1/619-641             UGAGUUCUCGAUCUCUAAAAUCG
        if not l.startswith('#') and not l.startswith('//'):
            seq = l.split()[1]
            print(seq)
