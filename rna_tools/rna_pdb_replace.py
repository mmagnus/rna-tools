#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Replace XYZ in file with insert.::

 ATOM     11  N1    A 2  27     303.441 273.472 301.457  1.00  0.00           N   # file
 ATOM      1  N1    A 2  27     300.402 273.627 303.188  1.00 99.99           N   # insert 
 ATOM     11  N1    A 2  27     300.402 273.627 303.188  1.00  0.00           N   # inserted
              XXXXXXXXXXXXXXXX    # part used to find lines to be replaced

 ATOM      1  P     A 2  27     295.653 270.783 300.135  1.00119.29           P   # next line

Example::

    python rna_pdb_replace.py input/t2-4-ACA.pdb input/to-replace.pdb

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')

from rna_tools.rna_tools_lib import RNAStructure, replace_atoms


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    parser.add_argument("insert", help="", default="") # nargs='+')    
    parser.add_argument("--output", help="by default <file>_rpl.pdb", default="") # nargs='+')    
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    f1 = args.file
    f2 = args.insert
    r1 = RNAStructure(f1)
    r2 = RNAStructure(f2)    

    t = replace_atoms(f1, f2, args.verbose)
    output = args.output
    if not output:
        output = args.file.replace('.pdb', '_rpl.pdb')

    with open(output, 'w') as f:
        f.write(t)
        print('Saved %s' % output)
        
    
            
