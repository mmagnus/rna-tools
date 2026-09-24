#!/usr/bin/env python
from __future__ import print_function
__docformat__ = 'reStructuredText'
import os
import Bio.PDB.PDBParser
import Bio.PDB.Superimposer
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import sqrt, array, asarray


BACKBONE_ATOMS = [
    'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"
]
# Atoms used to anchor base orientation per residue type
PURINE_ATOMS = ['N9', 'C8', 'C7']  # G, A
PYRIMIDINE_ATOMS = ['N1', 'C6', 'C2']  # C, U (and T)
PURINES = {'A', 'G', 'DA', 'DG'}
PYRIMIDINES = {'C', 'U', 'T', 'DT'}

ATOM_SELECTION_SUMMARY = "\n".join([
    f"BACKBONE_ATOMS = {BACKBONE_ATOMS}",
    "# Atoms used to anchor base orientation per residue type",
    f"PURINE_ATOMS = {PURINE_ATOMS}  # G, A",
    f"PYRIMIDINE_ATOMS = {PYRIMIDINE_ATOMS}  # C, U (and T)",
    f"PURINES = {sorted(PURINES)}",
    f"PYRIMIDINES = {sorted(PYRIMIDINES)}",
])


def get_atom_selection_summary():
    return ATOM_SELECTION_SUMMARY


def get_nucleotide_residues(struc):
    """Get nucleotides of the first model, in the file order.

    Waters, ions and ligands are skipped; modified nucleotides (HETATM with
    P and C1' atoms) are kept. The n-th residue returned here corresponds to
    the n-th residue of the sequence in an alignment (1-based)."""
    model = next(struc.get_models())
    residues = []
    for res in model.get_residues():
        hetflag = res.id[0]
        if hetflag == ' ' or (hetflag.startswith('H_') and res.has_id('P') and res.has_id("C1'")):
            residues.append(res)
    return residues


def get_sequence(struc):
    """Get the sequence of nucleotides (see get_nucleotide_residues), unknown residues as N."""
    seq = ''
    for res in get_nucleotide_residues(struc):
        resname = res.get_resname().strip().upper()
        if resname in ('A', 'C', 'G', 'U', 'T', 'DA', 'DC', 'DG', 'DT'):
            seq += resname[-1]
        else:
            seq += 'N'
    return seq

class RNAmodel:
    """RNAmodel

    :Example:

        >>> rna = RNAmodel("test_data/rp14/rp14_5ddp_bound_clean_ligand.pdb", [1], False, None)
        >>> rna.get_report()
        "File:  rp14_5ddp_bound_clean_ligand.pdb  # of atoms: 1 \\nresi:  1  atom:  <Atom C3'> \\n"

    :param fpath: file path, string
    :param residues: list of residues to use (and since we take only 1 atom, C3', this equals to number of atoms.
    :param save: boolean, save to output_dir or not
    :param output_dir: string, if save, save segments to this folder
    """
    #:returns: None
    #:rtype: None
    #"""
    def __init__(self, fpath, residues, save=False, output_dir=""):

        # parser 1-5 -> 1 2 3 4 5
        self.struc = RNAmodel.parse(fpath)
        self.fpath = fpath
        self.fn = os.path.basename(fpath)
        self.residues = residues #self.__parser_residues(residues)
        self.__get_atoms()
        #self.atoms = []
        if save:
            self.save(output_dir) # @save

    @staticmethod
    def parse(fpath):
        """Parse a PDB file, return Bio.PDB structure"""
        return Bio.PDB.PDBParser(QUIET=True).get_structure('', fpath)

    def __parser_residues(self, residues):
        """Get string and parse it
        '1 4 5 10-15' -> [1, 4, 5, 10, 11, 12, 13, 14, 15]"""
        rs = []
        for r in residues.split():
            l = parse_num_list(r)
            for i in l:
                if i in rs:
                    raise Exception('You have this resi already in your list! See', residues)
            rs.extend(l)
        return rs

    def __get_atoms(self):
        """Select atoms of the residues.

        self.residues are positions (1-based) in the sequence of the structure
        (see get_nucleotide_residues), not the residue numbers from the PDB file.
        Atoms are keyed by (index in self.residues, atom name) so two models
        built from lists of paired positions can be compared directly."""
        self.atoms = []
        self.atom_ids = []  # (index in self.residues, atomName)
        self.atom_lookup = {}
        nts = get_nucleotide_residues(self.struc)
        self.selected_residues = []
        for i, pos in enumerate(self.residues):
            if pos < 1 or pos > len(nts):
                raise Exception('residue position %s out of range (%s has %s nucleotides)' % (pos, self.fn, len(nts)))
            res = nts[pos - 1]
            self.selected_residues.append(res)
            # backbone atoms
            for atom_name in BACKBONE_ATOMS:
                self._append_atom_if_present(res, i, atom_name)
            # base atoms
            for atom_name in self._get_base_atoms(res):
                self._append_atom_if_present(res, i, atom_name)
        if len(self.atom_ids) <= 0:
            raise Exception('problem: no atoms were selected!: %s' % self.fn)
        return self.atoms

    def _append_atom_if_present(self, residue, index, atom_name):
        if residue.has_id(atom_name):
            atom = residue[atom_name]
            key = (index, atom_name)
            self.atom_lookup[key] = atom
            self.atom_ids.append(key)
            self.atoms.append(atom)
        else:
            pass  # silently ignore missing atoms

    def _get_base_atoms(self, residue):
        resname = residue.get_resname().strip().upper()
        if resname in PURINES or (resname and resname[0] in ('A', 'G')):
            return PURINE_ATOMS
        if resname in PYRIMIDINES or (resname and resname[0] in ('C', 'U', 'T')):
            return PYRIMIDINE_ATOMS
        return []

    def __str__(self):
        return self.fn #+ ' # beads' + str(len(self.residues))

    def __repr__(self):
        return self.fn #+ ' # beads' + str(len(self.residues))

    def get_report(self):
        """Str a short report about rna model"""
        t = ' '.join(['File: ', self.fn, ' # of atoms:', str(len(self.atoms)), '\n'])
        for r, a in zip(self.residues, self.atoms):
            t += ' '.join(['resi: ', str(r), ' atom: ', str(a), '\n'])
        return t

    def get_rmsd_to(self, other_rnamodel, output='', dont_move=False):
        """Calc rmsd P-atom based rmsd to other rna model"""
        sup = Bio.PDB.Superimposer()

        paired_self, paired_other = self._get_matched_atom_lists(other_rnamodel)
        if dont_move:
            coords = array([a.get_vector().get_array() for a in paired_self])
            other_coords = array([a.get_vector().get_array() for a in paired_other])
            s = SVDSuperimposer()
            s.set(coords, other_coords)
            return s.get_init_rms()

        sup.set_atoms(paired_self, paired_other)

        rms = round(sup.rms, 3)
        
        if output:
            io = Bio.PDB.PDBIO()
            sup.apply(self.struc.get_atoms())
            io.set_structure( self.struc )
            io.save("aligned.pdb")

            io = Bio.PDB.PDBIO()
            sup.apply(other_rnamodel.struc.get_atoms())
            io.set_structure( other_rnamodel.struc )
            io.save("aligned2.pdb")
        return rms

    def _get_matched_atom_lists(self, other_rnamodel):
        other_lookup = other_rnamodel.atom_lookup
        common_keys = [key for key in self.atom_ids if key in other_lookup]
        if not common_keys:
            raise Exception('No common atoms found between %s and %s' % (self.fn, other_rnamodel.fn))
        self_atoms = [self.atom_lookup[key] for key in common_keys]
        other_atoms = [other_lookup[key] for key in common_keys]
        return self_atoms, other_atoms

    def save(self, output_dir, verbose=True):
        """Save structures and motifs """
        folder_to_save =  output_dir + os.sep # ugly hack 'rp14/'
        try:
            os.makedirs(folder_to_save)
        except OSError:
            pass

        try:
            os.mkdir(folder_to_save + 'structures')
        except OSError:
            pass

        try:
            os.mkdir(folder_to_save + 'motifs')
        except OSError:
            pass

        RESI = set(id(r) for r in self.selected_residues)
        if not self.struc:
            raise Exception('self.struct was not defined! Can not save a pdb!')

        class BpSelect(Select):
            def accept_residue(self, residue):
                if id(residue) in RESI:
                    return 1
                else:
                    return 0

        io = PDBIO()
        io.set_structure(self.struc)
        fn = folder_to_save + 'structures' + os.sep + self.fn #+ '.pdb'
        io.save(fn)
        if verbose:
            print('    saved to struc: %s ' % fn)

        io = PDBIO()
        io.set_structure(self.struc)
        fn = folder_to_save +  'motifs/' + os.sep + self.fn #+ self.fn.replace('.pdb', '_motif.pdb')# #+ '.pdb'
        io.save(fn, BpSelect())
        if verbose:
            print('    saved to motifs: %s ' % fn)

#main
if __name__ == '__main__':
    import doctest
    doctest.testmod()

    #rna = RNAmodel("test_data/rp14/rp14_5ddp_bound_clean_ligand.pdb", [1], False, Non3e)
    #print(rna.get_report())
    a = RNAmodel("test_data/GGC.pdb", [46,47,48])
    b = RNAmodel("test_data/GUC.pdb", [31, 32, 33])

    print(a.get_rmsd_to(b))
    print(a.get_rmsd_to(b, dont_move=True))
    
