#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
import math
from numpy import log as ln

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-T", default = 310)    
    parser.add_argument("pkd", help="in M", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    from decimal import Decimal
    pkd = float(args.pkd)
    #print('pKd: %.2E' % Decimal(kd))
    kd = 10**-(pkd)
    print('Kd', kd)
    dG = args.T * 1.98720425864083*10**-3 * ln(kd)
    print('dG', round(dG, 2), 'kcal/mol')    
