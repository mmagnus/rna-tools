#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
import math

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("kd", help="in M", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    from decimal import Decimal
    kd = float(args.kd)
    print('Kd: %.2e' % Decimal(kd))
    print('pKd', -math.log10(kd))    
