#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import os

from pdb_parser_lib import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser('yapdb_parser')
    parser.add_argument('-r', '--report', help='get report',
                        action='store_true')
    parser.add_argument('--rosetta2generic', help='convert ROSETTA-like format to generic pdb',
                        action='store_true')
    parser.add_argument('-c', '--clean', help='get clean structure',
                        action='store_true')

    parser.add_argument('file', help='file') 
    parser.add_argument('outfile', help='outfile')   

    args = parser.parse_args()

    s = StrucFile(args.file)
    if args.report:
        print s.get_report()
        print s.get_preview()

    s = StrucFile(args.file)
    if args.clean:
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        print s.get_preview()
        s.write(args.outfile)

    if args.rosetta2generic:
        s = StrucFile(args.file)
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.fix_op_atoms()
        s.renum_atoms()
        print s.get_preview()
        s.write(args.outfile)
