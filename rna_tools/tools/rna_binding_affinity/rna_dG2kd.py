#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
calc Kd [M]
"""
from __future__ import print_function

import argparse
import math

from decimal import Decimal
from numpy import log as ln
from numpy import exp
from si_prefix import si_format


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-T", default = 310)
    parser.add_argument("dG", help="in kcal/mol", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    x = float(args.dG) / (args.T * 1.98720425864083*10**-3)
    kd = math.e ** x
    print('Kd:', kd, 'M')
    #print(si_format(kd, precision=2)) # 3.59 f
    print('Kd: %.2e' % Decimal(kd))
