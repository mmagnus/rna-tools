#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""rna_pdb_edit_occupancy_bfactor.py - a script to edit occupancy or bfactor in a PDB file.

Example::

   rna_pdb_edit_occupancy_bfactor.py --occupancy --select A:1-40,B:1-22 \\
                                    --set-to 0 \\
                                    19_Bujnicki_Human_4_rpr_n0-000001.pdb 
"""

from __future__ import print_function

from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.Atom import PDBConstructionWarning

import sys
import warnings
warnings.simplefilter('ignore', PDBConstructionWarning)

import re
import string
import argparse

def edit_occupancy_of_pdb(txt, pdb, pdb_out, bfactor, occupancy, set_to, set_not_selected_to, v=False):
    """Make all atoms 1 (flexi) and then set occupancy 0 for seletected atoms.

    Load the structure, and first set everything to be `set_not_selected_to`
    and then  set selected to `sel_to`.

    Args:

       txt (str): A:1-10, selection, what to change
       pdb (str): file to read as an input
       pdb_out (str): file to save an output
       bfactor (bool): if edit bfactor
       occupancy (bool): if edit occupancy
       set_to (float): set to this value, if within selection
       set_not_selected_to (float): set to this value, if not within selection
       v (bool): be verbose 

    Returns:

       bool: if OK, save an output to pdb_out

    """
    struc = PDB.PDBParser().get_structure('struc', pdb)

    # first loop to set everything to set_not_selected_to
    for s in struc:
         for c in s:
             for r in c:
                 for a in r:
                     if bfactor:
                         a.set_bfactor(set_not_selected_to)  
                     if occupancy:
                         a.set_occupancy(set_not_selected_to)

    # parser selected
    txt = txt.replace(' ','')
    if v: print(txt)
    l = re.split('[,:;]', txt)
    if v: print(l)

    # second loop to set_to for selected
    for i in l: # ['A', '1-10', '15', '25-30', 'B', '1-10']
        if i in string.ascii_letters:
            if v: print('chain', i)
            chain_curr = i
            continue

        if i.find('-') > -1:
            start, ends = i.split('-')
            if v: print(start, ends)
            if start > ends:
                print('Error: range start > end ' + i, file=sys.stderr)
                return False
            index = list(range(int(start), int(ends)+1))
        else:
            index=[int(i)]

        for i in index:
            if v: print(index)
            try:
                atoms = struc[0][chain_curr][i]
            except KeyError:
                if i == chain_curr:
                    print('Error: Chain ' + chain_curr + ' not found in the PDB structure', file=sys.stderr)
                else:
                    print('Error: Residue ' + chain_curr + ':' + str(i) + ' found in the PDB structure', file=sys.stderr)
                    return False

            for a in atoms:
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
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--bfactor', action='store_true', help="set bfactor")
    group.add_argument('--occupancy', action='store_true', help="set occupancy")

    parser.add_argument('--select', help='get chain, e.g A:1-10, works also for multiple chains, e.g A:1-40,B:1-22')
    parser.add_argument('--set-to', help='set value to', default=1)
    parser.add_argument('--set-not-selected-to', help='set value to', default=0)
    parser.add_argument('-o', '--output', help='file output')
    parser.add_argument('file', help='file')
    parser.add_argument('--verbose', action='store_true', help="be verbose")
    return parser

# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    if not args.output:
        args.output = args.file.replace('.pdb', '_out.pdb')
    edit_occupancy_of_pdb(args.select, args.file, args.output, args.bfactor, args.occupancy,
                              float(args.set_to), float(args.set_not_selected_to), args.verbose)
