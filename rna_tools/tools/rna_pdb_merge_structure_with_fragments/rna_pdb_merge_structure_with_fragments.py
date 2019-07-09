#!/usr/bin/env python

from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.Atom import PDBConstructionWarning

import warnings
warnings.simplefilter('ignore', PDBConstructionWarning)
import string
import re
import sys
import copy


class Struc:
    def __init__(self, ring, pdb, fragments, pdb_out, v):
        self.fragments = fragments
        self.ring = PDB.PDBParser().get_structure('ring', ring)
        self.pdb = PDB.PDBParser().get_structure('pdb', pdb)
        self.pdb_out = pdb_out
        self.v = v

    def detect_what_is_in_pdb(self):
        fragments = ''
        for c in self.pdb[0]:
            fragments += ';' + c.get_id() + ':'
            for r in c:
                fragments += str(r.get_id()[1]) + ','
            fragments = fragments[:-1]
        fragments = fragments[1:]  # to remoev first ;

        self.fragments = fragments
        print('fragments in pdb', self.fragments)

    def merge(self):
        v = self.v
        txt = self.fragments
        txt = txt.replace(' ', '')
        if v:
            print(txt)
        l = re.split('[,:;]', txt)
        if v:
            print(l)

        if v:
            print('ring', self.ring)

        for i in l:  # ['A', '1-10', '15', '25-30', 'B', '1-10']
            if i in string.ascii_letters:
                if v:
                    print('chain', i)
                chain_curr = i
                continue

            if i.find('-') > -1:
                start, ends = i.split('-')
                if start > ends:
                    return 'Error: range start > end ' + i
                index = list(range(int(start), int(ends) + 1))
            else:
                index = [int(i)]

            for i in index:
                try:
                    residue_ring = self.ring[0][chain_curr][i]
                except KeyError:
                    if i == chain_curr:
                        return 'Error: Chain ' + chain_curr + ' not found in the PDB structure (seq)'
                    else:
                        return 'Error: Residue ' + chain_curr + ':' + str(i) + ' found in the PDB structure (seq)'
                try:
                    residue_pdb = self.pdb[0][chain_curr][i]
                except KeyError:
                    if i == chain_curr:
                        return 'Error: Chain ' + chain_curr + ' not found in the PDB structure (seq)'
                    else:
                        return 'Error: Residue ' + chain_curr + ':' + str(i) + ' found in the PDB structure (seq)'

                if residue_ring.get_resname() == residue_pdb.get_resname():
                    if v:
                        print(i, residue_ring.get_resname(), residue_pdb.get_resname(), ' ...coping')
                    to_remove = self.ring[0][chain_curr][i].copy()
                    for a in to_remove:
                        self.ring[0][chain_curr][i].detach_child(a.get_id())  # remove atoms

                    for a in self.pdb[0][chain_curr][i]:
                        try:
                            self.ring[0][chain_curr][i].add(a)
                        except PDB.PDBExceptions.PDBConstructionException:
                            return 'Error: duplicated atoms in ' + chain_curr + ':' + str(i)
                else:
                    return 'Error: Residue ' + chain_curr + ':' + str(i) + ' not the same as in the sequence'

        io = PDBIO()
        io.set_structure(self.ring)
        io.save(self.pdb_out)
        print('Save as ', self.pdb_out)
        return False


if __name__ == '__main__':
    s = Struc(ring='test_data/1xjr_simrna_nstep_1.pdb',
              pdb='test_data/rna_sq_renumber.pdb', fragments='', pdb_out='tmp.pdb', v=True)
    s.detect_what_is_in_pdb()

    s = Struc(ring='test_data/1xjr_simrna_nstep_1.pdb', pdb='test_data/1xjr_root.pdb',
              fragments='A:1-14,37-45', pdb_out='test_data/tmp.pdb', v=True)
    err = s.detect_what_is_in_pdb()
    err = s.merge()
    if err:
        print(err)
