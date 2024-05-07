#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author: Guillaume Bouvier -- guillaume.bouvier@pasteur.fr
https://research.pasteur.fr/en/member/guillaume-bouvier/
2018-01-24 15:23:33 (UTC+0100)
https://bougui505.github.io/2018/01/24/compute_the_global_distance_test_(gdt)_with_pymol.html
"""

from pymol import cmd
import numpy as np

print("Loaded calc_gdt.py: gdt(sel1, sel2), gdt_hs(sel1, sel2)")
@cmd.extend
def gdt(sel1, sel2):
    """
    Calculate the GDT (Global Distance Test) between two RNA (using P atoms) selections.
    Args:
        sel1 (str): selection string for the protein to compare.
        sel2 (str): selection string for the reference protein.
    """
    # Select only C-alpha atoms
    sel1 += ' and name P'
    sel2 += ' and name P'
    cutoffs = [1.0, 2.0, 4.0, 8.0]
    # Align without transforming for raw alignment data
    cmd.align(sel1, sel2, cycles=0, transform=0, object='aln')
    # Get raw alignment data
    mapping = cmd.get_raw_alignment('aln')
    distances = []
    for pair in mapping:
        atom1 = '%s and id %d' % (pair[0][0], pair[0][1])
        atom2 = '%s and id %d' % (pair[1][0], pair[1][1])
        dist = cmd.get_distance(atom1, atom2)
        distances.append(dist)
    distances = np.array(distances)
    gdts = [(distances <= cutoff).mean() for cutoff in cutoffs]
    out = np.array(list(zip(cutoffs, gdts))).flatten()
    print("GDT_%d: %.4f; GDT_%d: %.4f; GDT_%d: %.4f; GDT_%d: %.4f;" % tuple(out))
    GDT = np.mean(gdts)
    print("GDT: %.4f" % GDT)
    # Color by distance
    cmd.spectrum('b', 'green_yellow_red', sel1)
    return GDT

@cmd.extend
def gdt_ha(sel1, sel2):
    """
    Calculate the GDT (Global Distance Test) High Accuracy between two RNA (using P atoms) selections.
    Args:
        sel1 (str): selection string for the protein to compare.
        sel2 (str): selection string for the reference protein.
    """
    # Select only C-alpha atoms
    sel1 += ' and name P'
    sel2 += ' and name P'
    cutoffs = [0.5, 1, 2, 4]
    # Align without transforming for raw alignment data
    cmd.align(sel1, sel2, cycles=0, transform=0, object='aln')
    # Get raw alignment data
    mapping = cmd.get_raw_alignment('aln')
    distances = []
    for pair in mapping:
        atom1 = '%s and id %d' % (pair[0][0], pair[0][1])
        atom2 = '%s and id %d' % (pair[1][0], pair[1][1])
        dist = cmd.get_distance(atom1, atom2)
        distances.append(dist)
    distances = np.array(distances)
    gdts = [(distances <= cutoff).mean() for cutoff in cutoffs]
    out = np.array(list(zip(cutoffs, gdts))).flatten()
    print("GDT_HS%d: %.4f; GDT_%d: %.4f; GDT_%d: %.4f; GDT_%d: %.4f;" % tuple(out))
    GDT = np.mean(gdts)
    print("GDT_HS: %.4f" % GDT)
    # Color by distance
    cmd.spectrum('b', 'green_yellow_red', sel1)
    return GDT

if __name__ == '__main__':
    # Assume the first two objects are the models to compare
    model1, reference = cmd.get_object_list('all')[:2]
    gdt_gs(model1, reference)
