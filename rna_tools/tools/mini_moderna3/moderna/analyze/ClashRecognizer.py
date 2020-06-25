#!/usr/bin/env python
#
# ClashRecognizer.py
#
# Class finding clashing residues.
#
# http://iimcb.genesilico.pl/moderna/
#

"""
One class is implemented here - ClashRecognizer.
It can detect clashes between different RNA residues
(not between atoms in the same residue).
"""

__author__ = "Tomasz Puton"
__copyright__ = "Copyright 2008, The Moderna Project"
__contributors__ = "Magdalena Rother, Kristian Rother"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

from Bio.PDB import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from rna_tools.tools.mini_moderna3.moderna.analyze.ChainConnectivity import are_residues_connected
from rna_tools.tools.mini_moderna3.moderna.Constants import ATOM_RADII, SEARCH_RADIUS

class ClashRecognizer:
    """
Class that can detect clashes between two different
RNA residues (not between atoms within a residue).

it implements three public methods:

- find_clashes_in_pdb(pdb_file)
- find_clashes_in_structure(struct)
- find_clashes_in_residues(residues)
    """    
    #TODO: atom radii and search radius as configurable parameters?
    
    radius_cache = {} # memorize radii by atom full names

    def compare_residues(self, one, two):
        # !!! KR: this should go into ModernaResidue
        return cmp(one.id[1], two.id[1])    

    def __get_model_atoms__(self, struct):
        """
        Generates lists of all atoms for each model in 
        a structure. 'struct' has to be a Structure object
        returned by PDBParser from Bio.PDB.PDBParser.
        """        
        for model in struct.get_list():
            list_of_atoms = []
            for atom in model.get_atoms():
                list_of_atoms.append(atom)
            yield list_of_atoms

    def __get_atom_radius__(self, atom):
        """
        Returns atom radius (float).
        'atom' has to be an Atom object returned by PDBParser
        """
        #TODO: make a dictionary for {atom name:radius} as a class 
        # variable to avoid multiple calculations.
        element = atom.fullname
        radius = self.radius_cache.get(element)
        if radius:
            return radius
        else:
            for char in element.strip()[:2]:
                if char in 'NCPOS':
                    radius = ATOM_RADII[char]
                    self.radius_cache[atom.fullname] = radius
                    return radius
            return 1.05
            # assume that the covalent radius for an unknown atom is 1.05 A
        #TODO: this could be a place for the memoize pattern.
    
    def __read_pdb_structure__(self, pdb_filename):
        """
        Makes a Structure object from a pdbfile
        """
        # KR: this probably can be outsourced to another module.
        parser = PDBParser()
        struct = parser.get_structure('',pdb_filename)
        return struct

    def __make_structure_from_residues__(self, residues):
        """
        Makes a Structure object either from a pdbfile or a list of residues
        """
        # KR: this probably can be outsourced to another module.
        struct = Structure('s')
        model = Model('m')
        n_chain = 1
        chain = Chain('c%i'%n_chain)

        for residue in residues:
            if chain.has_id(residue.id):
                model.add(chain)
                n_chain+=1
                chain = Chain('c%i'%n_chain)
            chain.add(residue)                   
            
        model.add(chain)
        struct.add(model)
        return struct

    def __find_neighbors__(self, struct):
        """Generates tuples of Atom objects that are neighbors."""
        model_atoms = self.__get_model_atoms__(struct)
        neighbors = [NeighborSearch(atom).search_all(SEARCH_RADIUS,'A') \
            for atom in model_atoms]
        return neighbors[0]
        #yield neighbors
        #KR: dropped the generator because it generated just one list of tuples.

    def __are_there_any_clashes__(self, neighbors):
        """
        The method eats output of the search_all method of NeighborSearch
        object. It generates tuples of clashing residues.
        """
        if len(neighbors) != 0:
            for first_atom, second_atom in neighbors:
                first_atom_radius = self.__get_atom_radius__(first_atom)
                second_atom_radius = self.__get_atom_radius__(second_atom)

                distance_between_atom_centers = first_atom - second_atom
               
                resi_1 = first_atom.get_parent()
                resi_2 = second_atom.get_parent()
                
                if (first_atom.get_name() == 'O3\''\
                and second_atom.get_name() == 'P') or \
                (first_atom.get_name() == 'P'\
                and second_atom.get_name() == 'O3\'') \
                and \
                (are_residues_connected(resi_1, resi_2) or \
                 are_residues_connected(resi_2, resi_1)):
                    # 15 May 2009:
                    # TP: feature requested by MM & KR
                    # Its really a temporary solution and it will be changed
                    # when ModernaStructure.are_residues_connected() is ready!
                    # (because we have to make sure that those two residues
                    # are linear neighbours)
                    minimum_distance_between_atom_centers = 1.413
                else:
                    minimum_distance_between_atom_centers = \
                    first_atom_radius + second_atom_radius
                
                #TODO: could be checked before calculating the atom-atom 
                # distance explicitly
                if resi_1 != resi_2 and distance_between_atom_centers < \
                    minimum_distance_between_atom_centers:
                    #We have a clash!
                    if self.compare_residues(resi_1, resi_2) <0:
                        yield resi_1, resi_2
                    else: 
                        yield resi_2, resi_1

    def __find_clashes__(self, struct):
        """
Finds clashes. Input can be PDB, Structure object and list.
Returns a list of clashing pairs of residues.
        """
        clashes = []
        neighbors = self.__find_neighbors__(struct)
        for clash in self.__are_there_any_clashes__(neighbors):
            if clash not in clashes:
                clashes.append(clash)
        return clashes
       
    def find_clashes_in_pdb(self, pdb_file):
        """
find_clashes_in_pdb(pdb_file)
- runs clash recognition function on a PDB file
- input is a path to a PDB file
- returns True if a clash is found and False if not
        """
        struct = self.__read_pdb_structure__(pdb_file)
        return self.__find_clashes__(struct)

    def find_clashes_in_structure(self, struct):
        """
find_clashes_in_structure(struct)
- runs clash recognition function on a Structure object returned by
  PDBParser from Bio.PDB.PDBParser
- returns True if a clash is found and False if not
        """
        return self.__find_clashes__(struct)
    
    def find_clashes_in_residues(self, residues):
        """
find_clashes_in_residues(self, residues)
- runs clash recognition function on a list of residues;
  a residue is an object of returned by PDBParser from 
  Bio.PDB.PDBParser
- returns True if a clash is found and False if not
        """
        if len(residues) == 0:
            return False
        struct = self.__make_structure_from_residues__(residues)
        return self.__find_clashes__(struct) 
