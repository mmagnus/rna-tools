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


def convert_cif_to_pdb(cif_file, add_header_to_output=True, version='', verbose=True,
                       split_conflicting_chains=True, output_path=None):
    """Convert an mmCIF file into one or more PDB files.

    Returns a list of generated PDB filenames. When a direct conversion fails
    (e.g. multi-character chain ids) the chains are split into separate models
    following the legacy CLI behavior.
    """
    from Bio.PDB import MMCIFParser, PDBIO

    def _prepend_header(pdb_file, remarks=None):
        if not add_header_to_output and not remarks:
            return
        new_content = ''
        if add_header_to_output:
            new_content += add_header(version) + '\n'
        if remarks:
            new_content += '\n'.join(remarks) + '\n'
        with open(pdb_file, 'r') as fh:
            new_content += fh.read()
        with open(pdb_file, 'w') as fh:
            fh.write(new_content)

    if output_path and split_conflicting_chains:
        raise ValueError('output_path is only supported when split_conflicting_chains is False')

    outputs = []
    parser = MMCIFParser()
    structure = parser.get_structure("structure_id", cif_file)
    pdb_file = output_path or cif_file.replace('.cif', '_fCIF.pdb')

    io = PDBIO()
    io.set_structure(structure)
    try:
        io.save(pdb_file)
        _prepend_header(pdb_file)
        outputs.append(pdb_file)
        if verbose:
            print(f'saved: {pdb_file}')
        return outputs
    except Exception:
        if verbose:
            print('Warning: some of the chains in this mmCIF file have multi-character ids.')

    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB import Structure, Model, PDBIO

    def has_high_rna_content(chain, threshold=0.8):
        rna_nucleotides = ['A', 'C', 'G', 'U', 'X']
        total_residues = 0
        rna_residues = 0
        for residue in chain:
            total_residues += 1
            if residue.get_resname().strip() in rna_nucleotides:
                rna_residues += 1
        if total_residues == 0:
            return False
        return (rna_residues / total_residues) >= threshold

    parser = MMCIFParser()
    structure = parser.get_structure("structure", cif_file)

    if not split_conflicting_chains:
        import string
        letters = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)
        letter_iter = iter(letters)
        for model in structure:
            for chain in model:
                try:
                    chain.id = next(letter_iter)
                except StopIteration:
                    raise RuntimeError('Too many chains to assign unique single-letter IDs in %s' % cif_file)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)
        _prepend_header(pdb_file)
        outputs.append(pdb_file)
        if verbose:
            print(f'saved: {pdb_file} (renamed chain ids)')
        return outputs

    import string
    letters = list(string.ascii_uppercase)

    for model in structure:
        for chain in model:
            if has_high_rna_content(chain):
                new_structure = Structure.Structure("new_structure")
                new_model = Model.Model(0)
                new_structure.add(new_model)

                chain_id_new = letters.pop(0)
                chain_id = chain.get_id()

                atom_count = 0
                for residue in chain:
                    for atom in residue:
                        atom_count += 1

                remarks = [f'REMARK rna chain {chain.id} -> {chain_id_new}']
                pdb_file = cif_file.replace('.cif', f'_{chain_id}_n{chain_id_new}_fCIF.pdb')
                if verbose:
                    print(f'rna chain {chain.id} -> {chain_id_new} {pdb_file} # of atoms: {atom_count}')

                chain.id = chain_id_new
                new_model.add(chain)

                io = PDBIO()
                io.set_structure(new_structure)
                io.save(pdb_file)
                _prepend_header(pdb_file, remarks)
                outputs.append(pdb_file)

    return outputs

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
        convert_cif_to_pdb(
            cif_file,
            add_header_to_output=not args.no_hr,
            version=version,
            verbose=True,
        )
