#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
To add missing atom, run::

  $ python rna_pdb_rnapuzzle_ready.py --fix_missing_atoms input/ACGU_no_bases.pdb > output/ACGU_no_bases_fixed.pdb
"""

import argparse
import os
import time

from pdb_parser_lib import *

version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
version = version[1].strip()

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--renumber_residues', help='renumber residues', action='store_true')
    parser.add_argument('--fix_missing_atoms', help='fix missing atoms', action='store_true')
    parser.add_argument('file', help='file')
    parser.add_argument('--no_hr', help='do not insert the header into files',
                        action='store_true')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    
    s = StrucFile(args.file)
    s.decap_gtp()
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.fix_op_atoms()
    s.renum_atoms()
    #print s.get_preview()
    #s.write(args.outfile)
    if not args.no_hr:
        add_header(version=version)
    s.get_rnapuzzle_ready(args.renumber_residues, args.fix_missing_atoms)
    print s.get_text()
