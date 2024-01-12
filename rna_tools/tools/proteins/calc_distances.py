#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    from Bio.PDB import PDBParser
    from Bio.PDB import Selection

    for f in args.file:
        print(f)
        # Replace 'pdb_id' with your PDB ID or file
        parser = PDBParser()
        pdb_id = ''
        structure = parser.get_structure(pdb_id, f)

        for model in structure:
            for chain in model:
                residues = list(chain)
                for i in range(len(residues) - 1):
                    current_residue = residues[i]
                    next_residue = residues[i + 1]

                    # Assuming standard naming, the carbonyl oxygen is named 'O'
                    atom = 'CA'
                    if atom in current_residue and atom in next_residue:
                        oxygen1 = current_residue[atom]
                        oxygen2 = next_residue[atom]
                        distance = oxygen1 - oxygen2
                        print(f"Distance between {atom} in residue {current_residue.get_id()} and {atom} in residue {next_residue.get_id()} is {distance:.2f} Å")
