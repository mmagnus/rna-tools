#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
import math
from si_prefix import si_format


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-T", default = 310)    
    parser.add_argument("kd", help="can be in format '3.6*10**-15' or '3.6*10^-15' [M]", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    from decimal import Decimal
    from numpy import log as ln
    if '^' in args.kd:
        args.kd = args.kd.replace('^', '**')
    dG = args.T * 1.98720425864083*10**-3 * ln(eval(args.kd)) # 3.6*10**-7)  
    
    print('dG', round(dG, 2), 'kcal/mol')
    #print(si_format(kd, precision=0))
