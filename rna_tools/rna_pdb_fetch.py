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

from rna_tools.rna_tools_lib import edit_pdb, add_header, get_version
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', help='', action='version', version=version)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('--no-hr', help='do not insert the header into files',
                        action='store_true')
    parser.add_argument("file", help="", default="", nargs='+')
    return parser, version


if __name__ == '__main__':
    parser, version = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for cif_file in args.file:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser()
        structure = parser.get_structure("structure_id", cif_file)
        pdb_file = cif_file.replace('.cif', '_fCIF.pdb')

        try:
            # Save to PDB format
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_file)

            print(f'saved: {pdb_file}')
            # open a file add remarks
            new_file = ''
            with open(pdb_file, 'r') as f:
                if not args.no_hr:
                    new_file += add_header(version) + '\n'
                new_file += f.read()

            with open(pdb_file, 'w') as f:
                f.write(new_file)

        except:
            print('Warning: some of the chains in this mmCIF file has chain names with more char than 1, e.g. AB, and the PDB format needs single-letter code, e.g. A.')
            def has_high_rna_content(chain, threshold=0.8):
                # RNA nucleotides: A, C, G, U, and X (you can modify as needed)
                rna_nucleotides = ['A', 'C', 'G', 'U', 'X']
                total_residues = 0
                rna_residues = 0

                # Count the total number of residues and RNA-like residues
                for residue in chain:
                    total_residues += 1
                    if residue.get_resname().strip() in rna_nucleotides:
                        rna_residues += 1

                # Calculate the proportion of RNA residues
                if total_residues == 0:
                    return False  # Avoid division by zero if chain has no residues

                rna_percentage = rna_residues / total_residues

                # Check if the percentage of RNA residues is greater than or equal to the threshold (80% by default)
                return rna_percentage >= threshold

            from Bio.PDB.MMCIFParser import MMCIFParser
            from Bio.PDB import MMCIFParser, Structure, Model, Chain

            # Initialize the parser
            parser = MMCIFParser()

            # Parse the structure
            structure = parser.get_structure("structure", cif_file)

            # Create a list of single-letter chain identifiers
            import string
            letters = list(string.ascii_uppercase)

            for model in structure:
                for chain in model:
                    if has_high_rna_content(chain):
                        # New structure
                        new_structure = Structure.Structure("new_structure")
                        new_model = Model.Model(0)  # Create a new model
                        new_structure.add(new_model)  # Add the new model to the new structure

                        chain_id_new = letters.pop(0)
                        chain_id = chain.get_id()

                        atom_count = 0
                        for residue in chain:
                              for atom in residue:
                                   atom_count += 1

                        remarks = []
                        remarks.append(f'REMARK rna chain {chain.id} -> {chain_id_new}')

                        pdb_file = cif_file.replace('.cif', f'_{chain_id}_n{chain_id_new}_fCIF.pdb')
                        print(f'rna chain {chain.id} -> {chain_id_new} {pdb_file} # of atoms: {atom_count}')

                        chain.id = chain_id_new
                        new_model.add(chain)

                        io = PDBIO()
                        io.set_structure(new_structure)

                        io.save(pdb_file)
                        # open a file add remarks
                        new_file = ''
                        with open(pdb_file, 'r') as f:
                            if not args.no_hr:
                                new_file += add_header(version) + '\n'
                            if remarks:
                                new_file += '\n'.join(remarks) + '\n'
                            new_file += f.read()

                        with open(pdb_file, 'w') as f:
                            f.write(new_file)
