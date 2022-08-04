#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Works both of M1 and Ubuntu.
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
    parser.add_argument("--sep", help="", default="\t")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sep = args.sep
    if list != type(args.file):
        args.file = [args.file]

    # print('#%i' % len(args.file))
    for index, f in enumerate(args.file):
        wrapper = ClashScore()
        try:
            result = wrapper.run(f, args.verbose)
        except:
            result = 'error'     
        #print(str(index + 1) + sep + str(round((index/len(args.file)) * 100)) + '%' + sep + f + sep + str(result))
        print(f + sep + str(result))
        wrapper.cleanup()
