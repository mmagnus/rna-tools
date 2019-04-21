#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Obtain a list of interaction in an RNA molecule where "Interaction" is purely distance based (defined by --cuf-off). Later, you can use it to calculate distance based INF, dINF ;-).

Example::

    [mm] rna_calc_inf$ git:(master)$ ./rna_calc_dinf.py  test_output/1Y26.pdb
    X    13   X   14          bp   G   C                  WW_cis   1
    X    13   X   83          bp   G   C                  WW_cis   1
    X    13   X   82          bp   U   C                  WW_cis   1
    X    14   X   15          bp   C   G                  WW_cis   1
    X    14   X   83          bp   G   G                  WW_cis   1
    X    14   X   81          bp   G   G                  WW_cis   1
    X    14   X   82          bp   U   G                  WW_cis   1

use clarna_compare.py::

    [mm] rna_calc_inf$ ./rna_calc_dinf.py  test_output/1Y26.pdb > 1Y26.pdb.outCR
    [mm] rna_calc_inf$ clarna_compare.py -iref 1Y26.pdb.outCR -ichk 1Y26.pdb.outCR
    1Y26.pdb.outCR                               1Y26.pdb.outCR      1.000      0.000      1.000      1.000      1.000      1.000      1.000      1.000

You can use ``-d`` to get a list of all interacting bases, something like::

    draw_dists([(13, 14),(13, 83),...(82, 83)])

so you can plot all interacting bases:

.. image:: ../../rna_tools/tools/rna_calc_inf/doc/1y26_dinf.png

Mind, that draw_dists works on C2 atoms, that might be different from atoms detected with the program (e.g. different base atom could be detected to make an interaction).

.. image:: ../../rna_tools/tools/rna_calc_inf/doc/1y26_dinf2.png

"""
from __future__ import print_function
import Bio.PDB
import argparse

v = False # verbose
ATOMS = {'G' : "N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split(),
         'A' : "N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split(),
         'U' : "N1 C2 O2 N3 C4 O4 C5 C6".split(),
         'C' : "N1 C2 O2 N3 C4 N4 C5 C6".split()}

class Pair:
    def __eq__(self, other):
        """Check if Pair is equal.
        """
        if (self.aid == other.aid and self.bid == other.bid and self.achain == other.achain and self.bchain == other.bchain) \
            or (self.aid == other.bid and self.bid == other.aid and self.achain == other.bchain and self.bchain == other.achain):
            return True


    def __init__(self, a, b, aname, bname, achain, bchain):
        """
        aid = A10
        """
        self.a = a
        self.b = b
        self.aname = aname
        self.bname = bname
        self.achain = achain  # A
        self.bchain = bchain  # A
        self.aid = achain + str(a)  # A10
        self.bid = bchain + str(b)  # A10

    def __repr__(self):
        # return ("%s    %i   %s   %i          bp %s %s                  WW_cis   1" % (self.achain, self.a, self.bchain, self.b, self.aname, self.bname))
        return ("%s    %i   %s   %i          bp G C                   WW_cis   1" % (self.achain, self.a, self.bchain, self.b))   #  , self.aname, self.bname))

    def get_indexes(self):
        """ get_dists([[62, 71], """
        return(self.a, self.b)

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file', help="a PDB file")

    parser.add_argument("-d", "--draw-dists",
                        action="store_true", help="")

    parser.add_argument("-c", "--cut-off",
                        action="store_true", help="", default=5)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


def infx(filename, draw_dists, cut_off):
    parser = Bio.PDB.PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.

    structures = parser.get_structure('struc', filename)
    structure = structures[0] # 'structures' may contain several proteins in this case only one.

    atoms = Bio.PDB.Selection.unfold_entities(structure, 'A') # Atom
    ns = Bio.PDB.NeighborSearch(atoms)

    pairs = []

    # I'm working on A and B residues.
    for atom in atoms:
        if v: print(atom.get_parent().get_resname().strip(), atom.get_name())
        if atom.get_name() in ATOMS[atom.get_parent().get_resname().strip()]:  # '  C' so strip it
            if v: print('    Get distance')
            #target_atom = structure['A'][311]['CA']
            # Bs  = Bio.PDB.Selection.unfold_entities(structure, 'A')
            # print(atom.get_coord(), atom.get_parent())
            close_atoms = ns.search(atom.get_coord(), cut_off)
            for close_atom in close_atoms:

                bres = close_atom.get_parent()
                ares = atom.get_parent()

                if bres != ares:
                    a = bres.get_id()[1]  # to get res index
                    b = ares.get_id()[1]
                    # create an object Pair
                    p = Pair(b, a, bres.get_resname(), ares.get_resname(), bres.get_parent().id, ares.get_parent().id)
                    # and check if add it to the list?
                    if p not in pairs:
                        pairs.append(p)

        #print("A    %i   A   %i          bp C G                  X   1" % p)

    #print(close_atoms)
    #for a in close_atoms:
    #    print(a, a.get_parent())
    if not draw_dists:
        print('Classifier: Clarna')
        for p in pairs:
            print(p)

    if draw_dists:
        print('draw_dists([' + ','.join([str(p.get_indexes()) for p in pairs]) + '])')


# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    infx(args.file, args.draw_dists, args.cut_off)
