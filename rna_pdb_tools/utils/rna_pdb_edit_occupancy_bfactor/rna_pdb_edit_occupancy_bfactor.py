#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""rna_pdb_edit_occupancy_bfactor.py - edit occupancy or bfactor in PDB file.

Example::

   rna_pdb_edit_occupancy_bfactor.py --occupancy --select A:1-40,B:1-22 \
                                    --set-to 0 \
                                    19_Bujnicki_Human_4_rpr_n0-000001.pdb


   rna_pdb_edit_occupancy_bfactor.py --occupancy \
                                     --select A:1-2 \
                                     --select-atoms P+C4\' \
                                     --set-to 10 \
                                     -o test_data/3w3s_homologymodel_out.PD
                                     --set-not-selected-to 8
                                     test_data/3w3s_homologymodel.pdb

"""

from __future__ import print_function
import re
import string
import argparse

from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.Atom import PDBConstructionWarning

import sys
import warnings
warnings.simplefilter('ignore', PDBConstructionWarning)


def edit_occupancy_of_pdb(txt, pdb, pdb_out, bfactor, occupancy, set_to,
                          set_not_selected_to, select_atoms, v=False):
    """Change ouccupancy or bfactor of pdb file.

    Load the structure, and first set everything to be `set_not_selected_to`
    and then  set selected to `sel_to`.

    Args:

       txt (str): A:1-10, selection, what to change
       pdb (str): filename to read as an input
       pdb_out (str): filename to save an output
       bfactor (bool): if edit bfactor
       occupancy (bool): if edit occupancy
       set_to (float): set to this value, if within selection
       set_not_selected_to (float): set to this value, if not within selection
       select_atoms (str): P, P+C4\\', use + as a separator
       v (bool): be verbose

    Returns:

       bool: if OK, save an output to pdb_out

    .. warning:: this function requires BioPython
    """
    struc = PDB.PDBParser().get_structure('struc', pdb)

    # first loop to set everything to set_not_selected_to
    if set_not_selected_to:
        for s in struc:
            for c in s:
                for r in c:
                    for a in r:
                        if bfactor:
                            a.set_bfactor(set_not_selected_to)
                        if occupancy:
                            a.set_occupancy(set_not_selected_to)

    # parser selected
    txt = txt.replace(' ', '')
    l = re.split('[,:;]', txt)
    if v:
        print(txt)
        print(l)

    # second loop to set_to for selected
    for i in l:  # ['A', '1-10', '15', '25-30', 'B', '1-10']
        if i in string.ascii_letters:
            if v:
                print('chain', i)
            chain_curr = i
            continue

        if '-' in i:
            start, ends = i.split('-')
            if v:
                print(start, ends)

            if int(start) > int(ends):
                print('Error: range start > end ' + i, file=sys.stderr)
                return False
            index = list(range(int(start), int(ends)+1))
        else:
            index = [int(i)]

        for i in index:
            if v:
                print(index)
            try:
                atoms = struc[0][chain_curr][i]
            except KeyError:
                if i == chain_curr:
                    print('Error: Chain ' + chain_curr +
                          ' not found in the PDB structure', file=sys.stderr)
                else:
                    print('Error: Residue ' + chain_curr + ':' + str(i) +
                          ' found in the PDB structure', file=sys.stderr)
                    return False

            for a in atoms:
                # if run it for all atoms
                if not select_atoms:
                    if bfactor:
                        a.set_bfactor(set_to)
                    if occupancy:
                        a.set_occupancy(set_to)
                else:  # if select_atoms
                    if '+' in select_atoms:
                            # only for striping atom names
                            select_atoms_list = [x.strip() for x in
                                                 select_atoms.split('+')]
                    else:
                            select_atoms_list = select_atoms
                    if a.id in select_atoms_list:
                        if bfactor:
                            a.set_bfactor(set_to)
                        if occupancy:
                            a.set_occupancy(set_to)

    io = PDBIO()
    io.set_structure(struc)
    io.save(pdb_out)
    print('Saved ', pdb_out)
    return True


def get_parser():
    """Get parser."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--bfactor', action='store_true', help="set bfactor")
    group.add_argument('--occupancy', action='store_true', help="set occupancy")

    p.add_argument('--select', help='get chain, e.g A:1-10, works also for multiple chains'
                                    'e.g A:1-40,B:1-22')
    p.add_argument('--set-to', help='set value to, default is 1', default=1)
    p.add_argument('--set-not-selected-to', help='set value to, default is 0', default=0)
    p.add_argument('-o', '--output', help='file output')
    p.add_argument('--verbose', action='store_true', help="be verbose")
    p.add_argument('--select-atoms', help="select only given atoms"
                                          "can be only one atom, e.g. P"
                                          "or more, use \\' for prims, e.g. P+C4\\'")
    p.add_argument('file', help='file')
    return p


# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    if not args.output:
        args.output = args.file.replace('.pdb', '_out.pdb')
    edit_occupancy_of_pdb(args.select, args.file, args.output, args.bfactor, args.occupancy,
                          float(args.set_to), float(args.set_not_selected_to),
                          args.select_atoms, args.verbose)
