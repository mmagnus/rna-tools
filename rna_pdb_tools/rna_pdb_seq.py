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

def sort_nicely( l ):
   """ Sort the given list in the way that humans expect.
   http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
   """
   convert = lambda text: int(text) if text.isdigit() else text
   alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
   l.sort( key=alphanum_key )
   return l

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('files', help='files', nargs='+')
    parser.add_argument('-v', '--verbose', help='',
                        action='store_true')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    for f in sort_nicely(args.files):
        s = StrucFile(f)
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        #print s.get_preview()
        if args.verbose:
            add_header(version=version)
        print s.get_seq()
        print
