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

from Bio.PDB import PDBParser, Superimposer
import numpy as np

# Function to parse a structure from a PDB file
def parse_structure(file_path):
    parser = PDBParser()
    return parser.get_structure('protein', file_path)

# Function to calculate GDT
def calculate_gdt(structure1, structure2, thresholds):
    # Superimpose structures
    ref_atoms = list(structure1.get_atoms())
    sample_atoms = list(structure2.get_atoms())
    
    # Ensuring both structures have the same number of atoms
    if len(ref_atoms) != len(sample_atoms):
        raise ValueError("Structures do not have the same number of atoms")
    
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_atoms)
    
    # Calculate distances after superimposition
    distances = np.array([ref_atom - sample_atom for ref_atom, sample_atom in zip(ref_atoms, sample_atoms)])
    
    # Calculate GDT for given thresholds
    gdt_scores = {}
    for threshold in thresholds:
        matches = np.sum(distances < threshold)
        gdt_scores[threshold] = (matches / len(ref_atoms)) * 100  # GDT percentage
    
    return gdt_scores


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-ha', "--high-accuracy", help="", default="", action="store_true")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("struc1", help="", default="") # nargs='+')
    parser.add_argument("struc2", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    # Example usage
    structure1 = parse_structure(args.struc1)
    structure2 = parse_structure(args.struc2)
    thresholds = [1, 2, 4, 8]  # Thresholds in Angstroms
    if args.high_accuracy:
        thresholds = [0.5, 1, 2, 4]
    gdt_scores = calculate_gdt(structure1, structure2, thresholds)
    print("GDT Scores:", gdt_scores)
