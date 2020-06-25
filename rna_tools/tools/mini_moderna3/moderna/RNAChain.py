#!/usr/bin/env python
#
# RNAChain.py
#
# Superclass that supports work with PDB.Structure objects.
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


from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from rna_tools.tools.mini_moderna3.moderna.RNAResidue import RNAResidue
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaAlphabet import alphabet
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence
from rna_tools.tools.mini_moderna3.moderna.Constants import MISSING_RESIDUE
from rna_tools.tools.mini_moderna3.moderna.util.Errors import RNAChainError, ModernaStructureError
from rna_tools.tools.mini_moderna3.moderna.analyze.ChainConnectivity import are_residues_connected
import os


class RNAChain:
    """
    Deals with chains from PDB.Structure objects. RNAChain allows to contain only one chain
    Is iterable as a list of all residues from this structure as ModernaResidue objects.
    Per m[x] residues can from structure m can be accessed by their numbers x.
    """

    def __init__(self, data_type=None, data=None, chain_name=None, seq=None, new_atoms=True):
        """
        Arguments:
        - data type - '' (to create a structure with empty chain), 'file', 'structure', or 'chain'
        - data - a pdb file name, a Bio.PDB.Structure instanse, a Bio.PDB.Chain instance
        - chain name
        - optional Sequence object, if the sequence is known (saves time)
        """
        # KR: I see 3 things which we could improve.
        #  - recognize the data type automatically
        #  - ModernaStructure is really a chain. Should this be changed?
        #  - give the ModernaStructure(PDB.Structure/Chain) another try.
        self.moderna_residues = {}
        self.chain_name = chain_name or 'A'

        if data_type == 'file':
            data = self._get_struct_from_file(data)
            
        if data_type in ['file', 'structure']:
            data = self._get_chain_from_struct(data, self.chain_name)
            
        if data_type in ['file', 'structure', 'chain', 'residues']:
            self._create_moderna_residues(data, seq, new_atoms)
        elif data_type: 
            raise RNAChainError('Unproper mode type.')
        # empty structures are OK.

        self.remove_empty_residues()    
        self.index = 0
        self.sort_residues()

    # 
    # methods helping with initializing object
    # 
    def _get_struct_from_file(self, filename):
        """Returns a Bio.PDB.Structure object."""
        if not os.access(filename, os.F_OK):
            raise RNAChainError('Structure file does not exist: %s'%filename)
        return PDBParser().get_structure('', filename)

    def _get_chain_from_struct(self, struct, chain_name):
        if chain_name not in struct[0].child_dict:
            raise RNAChainError("Chain '%s' does not exist."%chain_name)
        return struct[0][chain_name]

    def _create_moderna_residues(self, data, seq, new_atoms):
        """Creates dictionary of RNA residues from iterable of Bio.PDB.Residues."""
        if seq:
            for resi, aentry in zip(data, seq):
                temp = RNAResidue(resi, aentry, new_atoms=new_atoms)
                self.moderna_residues[temp.identifier] = temp
        else:
            for resi in data:
                temp = RNAResidue(resi, new_atoms=new_atoms)
                self.moderna_residues[temp.identifier] = temp
    # 
    # data management
    #
    def __bool__(self):
        """RNAChains are always True."""
        return True

    def __len__(self):
        """Returns number of residues."""
        return len(list(self.moderna_residues.keys()))
        
    def __repr__(self):
        return "<RNA structure; chain '%s'; %i residues>" % (self.chain_name, len(self))

    def __getitem__(self, args):
        """
        Allows to gain a residue (as a ModernaResidue instance) with given identifier
        or a dict of residues (as ModernaResidues instances) with given numbers (key - number, value - ModernaResidue) 
        """
        if type(args) == str:
            if args in list(self.moderna_residues.keys()): 
                return self.moderna_residues[args]
            else: 
                raise RNAChainError('There is no such residue: %s' %args)
        elif type(args) == slice:
            return self._get_residues_in_region(args.start, args.stop)
        else:
            raise RNAChainError('Bad argument type: %s'%str(args)) 

    def _get_residues_in_region(self, start_id=None,  stop_id=None):
        result = []
        if not start_id: 
            start_id = self.first_resi.identifier
        if not stop_id: 
            stop_id = self.last_resi.identifier
        self[start_id], self[stop_id] # check if residues are there
        all_resi = [resi for resi in self]
        resi_include = False
        for resi in all_resi:
            if resi.identifier == start_id: 
                resi_include = True
            if resi_include: 
                result.append(resi)
            if resi.identifier == stop_id: 
                resi_include = False
        return result
        
    def get_region(self, start_id=None,  stop_id=None):
        """Returns residues in the given range as an RNAChain."""
        subst = self._get_residues_in_region(start_id, stop_id)
        return self.__class__('residues', subst, self.chain_name)
    

    def __iter__(self):
        residues_list = []
        sorted_numbers = self.sort_residues()
        for num in sorted_numbers:
            residues_list.append(self[num])
        return residues_list.__iter__()
    

    def remove_residue(self, residue_identifier):
        """Removes the residue with the given identifier."""
        try:
            del self.moderna_residues[residue_identifier]
        except KeyError:
            raise ModernaStructureError('There is no residue with number %s. Could not remove residue.' % residue_identifier)
    
    def remove_empty_residues(self):
        for resi in self:
            if len(resi.child_list) == 0:
                self.remove_residue(resi.identifier)

    def get_all_atoms(self):
        """
        Returns a list of PDB.Atom objects, representing all atoms in the structure.
        """
        all_atoms = []
        for resi in self:
            all_atoms += [atom for atom in resi]
        return all_atoms   

    def get_structure(self, name='RNA chain'):
        """Returns chain as a PDB.Structure object."""
        struc = Structure(name)
        model = Model(0)
        chain = Chain(self.chain_name)
        struc.add(model)
        struc[0].add(chain)

        for resi in self:
            struc[0][self.chain_name].add(resi)
        return struc
  
  
    def write_pdb_file(self, output_file_name='structure.pdb'):
        """
        Writes the structure to a pdb file.
                
        Arguments:
        * output file name(/path) (by default 'structure.pdb')
        """
        struc = self.get_structure()
        output = PDBIO()
        output.set_structure(struc)
        output.save(output_file_name)

    def cmp_for_moderna_residues(self, x, y):
        """x, y have to be ModernaResidue objects!!!"""
        #TODO: move to RNAResidue
        x_nr = x.get_full_id()[0][1]
        x_ins_id = x.get_full_id()[0][2]
        y_nr = y.get_full_id()[0][1]
        y_ins_id = y.get_full_id()[0][2]
        
        if x_nr < y_nr:
            return -1
        if x_nr > y_nr:
            return 1
        if x_nr == y_nr:
            if x_ins_id < y_ins_id:
                return -1
            if x_ins_id > y_ins_id:
                return 1
            if x_ins_id == y_ins_id:
                return 0
     
    def sort_residues(self):
        """
        Sorts residues in a chain object (also negative numbers),
        so that they can be written in the right order.
        """
        from functools import cmp_to_key
        sorted_residues = sorted(list(self.moderna_residues.values()), key=cmp_to_key(self.cmp_for_moderna_residues))
        sorted_residue_numbers = [resi.identifier for resi in sorted_residues]
        return sorted_residue_numbers

    def get_sequence(self):
        """
        Returns a Sequence object containing the one-letter sequence
        of this ModernaStructure instance, including modifications.
        """
        seq = []
        previous_resi = None                

        for res in self:
            if previous_resi:
                if not are_residues_connected(previous_resi, res):
                    seq.append(alphabet[MISSING_RESIDUE])
            seq.append(alphabet[res.long_abbrev])
            previous_resi = res
        return Sequence(seq)

    @property
    def first_resi(self):
        """Returns the first residue."""
        return list(self)[0]

    @property
    def last_resi(self):
        """Returns the last residue."""
        return list(self)[-1]
