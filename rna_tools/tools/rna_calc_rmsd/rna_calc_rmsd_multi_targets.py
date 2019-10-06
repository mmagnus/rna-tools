#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

# run for one do it for all
# merge tables and calculcate statiics


"""
import argparse
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # parser.add_argument('-args', help="", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--models", help="", nargs='+')
    parser.add_argument("--targets", help="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    
    print(args.models)
    print(args.targets)

    for target in args.targets:
        cmd = "rna_calc_rmsd.py -t " + target + " " + ' '.join(args.models)
        os.system(cmd)

