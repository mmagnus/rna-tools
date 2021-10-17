#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import os
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    xx = "cHW tHW cHH tHH	cHS	tHS	cSW	tSW	cSH	tSH	cSS	tSS tsS".split()
    yy = "cWW tWW cHW tHW cSW tSW cWH tWH cHH tHH cSH tSH".split()
    c = 1
    for x in xx:
        for y in yy:
            cmd = 'wget http://rna.bgsu.edu/triples/script/triple_' + y + '_' + x + '.html'
            print(cmd)
            print(c)
            c += 1
            os.system(cmd)
