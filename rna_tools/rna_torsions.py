#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Examples::

    $ rna_torsions.py ./input/4GXY_min.pdb
    f, alphaprime, beta
    input ./input/4GXY_min.pdb <Residue G het=  resseq=2 icode= >, -64.20924484900823, -143.18546007904766
    input ./input/4GXY_min.pdb <Residue C het=  resseq=3 icode= >, 2.3394112025736815, 70.4052871669199

Comparison::

    $ rna_x3dna.py input/4GXY_min.pdb -s
    input: input/4GXY_min.pdb
       nt id   res  alpha   beta  gamma  delta  epsilon  zeta      e-z           chi     phase-angle sugar-type  ssZp    Dp  splay     paired
    0   1  G  A.G2    NaN -143.2  153.7   82.5    -92.3 -31.9  -60(..)  -179.0(anti)  19.5(C3'-endo)  ~C3'-endo  4.39  4.56  18.32  no paired
    1   2  C  A.C3 -111.4   70.4  160.0   80.6      NaN   NaN      NaN  -177.6(anti)  11.1(C3'-endo)  ~C3'-endo   NaN   NaN    NaN  no paired

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import math
try:
    from Bio import PDB
except ImportError:
    print("Biopython is not detected. It is required for some functions.")

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]


    print(f'f, alphaprime, beta')    
    for f in args.file:
        #print(f'input {f}, ')

        parser = PDB.PDBParser()
        struct = parser.get_structure('', f)
        
        for m in struct:
              for c in m:
                    for r in c:
                        #print(r)#, end=' ')
                        # check what residues is this one       
                        x1, x2, x3, x4 = r["C5'"].get_vector(), r["O5'"].get_vector(), r['P'].get_vector(), r['OP1'].get_vector()
                        a = PDB.vectors.calc_dihedral(x4, x3, x2, x1) # !!! order
                        alphaprime = math.degrees(a)
                        x1, x2, x3, x4 = r["C5'"].get_vector(), r["O5'"].get_vector(), r['P'].get_vector(), r["C4'"].get_vector()
                        b = PDB.vectors.calc_dihedral(x3, x2, x1, x4)
                        beta = math.degrees(b) # % 360
                        print(f'input {f} {r}, {alphaprime}, {beta}')
