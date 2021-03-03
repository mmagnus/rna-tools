#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Standalone tool to run the RNAkb class in the terminal.

All atom mode does not really work, see the documentation of the RNAkb class.
"""
from __future__ import print_function

import argparse
import os
from rna_tools.tools.mq.RNAkb.RNAkb import RNAkb


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-p", "--potential", default="5pt", help="5pt")
    parser.add_argument("file", help="a PDB file, one or more", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

#    print('fn,Bond" "Angle" "Proper Dih." "Improper Dih." "LJ-14" "Coulomb-14" "LJ (SR)" "Coulomb (SR)" "Potential" "Kinetic En." "Total Energy" "Temperature" "Pressure (bar)"'.split()))

    print('fn,Bond,Angle,Proper Dih.,Improper Dih.,LJ-14,Coulomb-14,LJ (SR),Coulomb (SR),Potential,Kinetic En.,Total Energy,Temperature,Pressure (bar)')
    for f in args.file:
        wrapper = RNAkb(sandbox=True)
        result = wrapper.run(f, args.potential, args.verbose)
        print(','.join(result))
        #print(f + ',' + ','.join([str(x) for x in result]))
        #wrapper.cleanup()
        
