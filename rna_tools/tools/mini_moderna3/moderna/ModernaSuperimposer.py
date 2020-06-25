#!/usr/bin/env python
#
# ModernaSuperimposer.py
#
# Wrapper for Bio.PDB.Superimposer.
# 
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"


from Bio.PDB  import Superimposer
from rna_tools.tools.mini_moderna3.moderna.util.Errors import ModernaSuperimposerError


class ModernaSuperimposer:
    """
    Wrapper for Bio.PDB.Superimposer, that can handle
    extracting atoms with predefined names.
    Allows to rotate and translate atoms in space.
    """
    def __init__(self, fixed=None, moved=None, moved_atoms=None):
        """
        fixed - list of atoms the moved atom will be superimposed upon.
        moved - list of atoms that will be superimposed on the fixed ones.
        moved_atoms - all atoms whose coordinates should be transformed.
        If all three parameters are given, the superpositions 
        will be carried out immediately.
        """
        self.fixed = fixed
        self.moved = moved
        self.moved_atoms = moved_atoms
        self.rmsd = None
        if fixed and moved and moved_atoms: 
            self.superimpose()


    def superimpose(self):
        """
        Performs the superimposition.
        Returns RMSD. 
        """
        if not self.fixed or not self.moved: 
            raise ModernaSuperimposerError('There are no fixed or moved atoms. Can not calculate rotation and translation matrix.')
        if not self.moved_atoms: 
            raise ModernaSuperimposerError('There are no atoms for superimposition given. Can not applay rotation and translation matrix')
        sup = Superimposer()
        sup.set_atoms(self.fixed, self.moved)
        sup.apply(self.moved_atoms)
        self.rmsd = sup.rms
        return self.rmsd


    def get_atoms(self, resi_list, atom_names_list, mode='fixed'): 
        """
        Retrieves all atoms with given names from the given list of residues.
        Returns a list of PDB.Atom objects.
        Sets superimposition atoms - fihed or moved depending on given mode.

        Arguments:
        - resi_list - list of PDB.Residue objects
        - atom_names_list - list of atom names
        - mode - 'fixed' or 'moved' depending on which superposition atoms will be given
        """
        atoms_list = []
        for resi in resi_list:
            for name in atom_names_list:
                try:
                    atoms_list.append(resi[name])
                except KeyError:
                    raise ModernaSuperimposerError("Atom %s not available for superposition."%name)
                
        if mode == 'fixed': 
            self.fixed = atoms_list
        elif mode == 'moved': 
            self.moved = atoms_list
        return atoms_list

    def set_moved_atoms(self, atoms_list):
        self.moved_atoms = atoms_list
        #TODO: remove this method
