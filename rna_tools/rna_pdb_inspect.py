#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
from rna_tools_lib import RNAStructure
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
from Bio.PDB import PDBParser


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    res_std= {
        'G': "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split(), # 23 
        'A': "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split(),    # 22
        'U': "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split(),          # 20
        'C': "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split() }         # 20

    def is_rna_chain(chain):
        AMINOACID_CODES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                   "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
                   "TRP", "TYR", "VAL"]
        RES = ['A', 'G', 'U', 'C']

        aa = []
        na = []
        for r in chain:
            r = r.get_resname()
            if r in AMINOACID_CODES:
                aa.append(r)
            if r in RES:
                na.append(r)

        aa = float(len(aa)) / len(chain)
        na = float(len(na)) / len(chain)
        #ic([aa, na, len(chain)])
        if aa == 0 and na == 0:
            return 'error'
        if na > aa:
            return True
        return False

    n = 0
    for f in args.file:
        print(f'Input: {f}')
        r = RNAStructure(f)
        
        print(f'{len(r.get_all_chain_ids())} of chains:', r.get_all_chain_ids())
        # get_residue() # even as lines

        parser = PDBParser()
        struc = parser.get_structure('', f)
        for model in struc:
            for chain in model:
                print(f'  {chain.id} of RNA? {is_rna_chain(chain)}')
                if is_rna_chain(chain):
                    # is it protein?
                    for res in chain:
                        #print(f'    {res.get_resname()} {res}')
                        rtype = res.get_resname()
                        for i, a in enumerate(res):
                            try:
                                res_std[rtype][i]
                            except:
                                #print(f'     modification {rtype} or index out of range')
                                pass
                            else:
                                if res_std[rtype][i] != a.id:
                                    print(f'      !wrong order of atoms in residue ! {res_std[rtype][i]} {a.id}')
