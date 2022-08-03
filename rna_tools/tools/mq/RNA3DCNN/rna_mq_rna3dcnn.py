#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Works on Ubuntu, on M1 a problem with libs.
"""
from __future__ import print_function
import argparse
import os
from rna_tools.tools.mq.RNA3DCNN.RNA3DCNN import RNA3DCNN


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('-d', "--done", help="a csv with already done scores", default="")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]
    print('id,fn,rna3dcnn')
    for i, f in enumerate(args.file):
        wrapper = RNA3DCNN()
        result = wrapper.run(f, args.verbose)
        print(','.join([str(i + 1), os.path.basename(f), str(result)]))
