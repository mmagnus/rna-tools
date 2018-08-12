#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
from Seq import RNASequence
import argparse
import sys

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--cyclefold', action="store_true")

    parser.add_argument("seq", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if args.cyclefold:
        seq = RNASequence(args.seq)
        print(seq)
        print(seq.predict_ss('rnastructure_CycleFold')[1], seq.predict_ss('rnastructure_CycleFold')[0])
