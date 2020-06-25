#!/usr/bin/env python
#
# PdbConverter.py
#
# Supports conversion between old and new pdb files
# - ribose descriptor (* <---> ')
# - phosphate group nomenclature (O1P, O2P, O3P <---> OP1, OP2, OP3)
# - residue name position (18 <---> 19 <---> 20)
# - residue flag (ATOM <---> HETATM <---> ATOM/HETATM)
# 
# Converion can be carried out on:
# http://iimcb.genesilico.pl/moderna/ 
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"


import re,os

from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBIO
from .moderna import ModernaResidue


class ConversionParameters():
    def __init__(self,  ribose=None, phosphate=None, res_name=None, flag=None):
        self.ribose = ribose
        self.phosphate = phosphate
        self.res_name = res_name
        self.flag = flag
    

# superclass
class Converter():
    """
Superclass for ModernaStructureConverter and PdbStructureConverter.
    
Supports conversion between old and new pdb files:
- ribose descriptor (* <---> ')
- phosphate group nomenclature (O1P, O2P, O3P <---> OP1, OP2, OP3)
- residue name position (18 <---> 19 <---> 20)
- residue flag (ATOM <---> HETATM <---> ATOM/HETATM) 
    """

    def __init__(self, parameters):
        self.parameters = parameters
    
    
    def get_full_at_name(self, new_name):
        """
        Formats atom fullname.
        """
        new_name = new_name.strip()
        if len(new_name)<4: return ' '+new_name
        else: return new_name


    def change_at_name(self, at, new_name):
        """
        Changes single atom name.
        Returns a new Atom object.
        """
        at_fullname = self.get_full_at_name(new_name)
        element = new_name.strip()[0]
        if element not in "CHONSP":
            element = 'C'
        new_at = Atom(new_name, at.coord, at.bfactor, at.occupancy, at.altloc, at_fullname, at.serial_number, element=element)
        return new_at
        
        
    def get_new_id(self, resi):
        """
        Returns a new id with proper flag for given residue.
        (Flag means ATOM or HETATM)
        """
        old_id = resi.id
        # mark all resi as HETATM
        if self.parameters.flag == 'H': return ('H_'+resi.resname.strip(), old_id[1], old_id[2])
        # mark all resi as ATOM
        elif self.parameters.flag == 'A': return (' ', old_id[1], old_id[2])
        # marke standard resi as ATOMS and modifications as HETATM
        elif self.parameters.flag == 'HA':
            r = ModernaResidue(resi)
            if r.category == 'modified':  return ('H_'+resi.resname.strip(), old_id[1], old_id[2])
            elif r.category =='standard RNA': return (' ', old_id[1], old_id[2])
            else: return old_id
        else: return old_id
        
        
    def get_new_name(self, resi):
        """
        Returns a new name for given residue.
        (Base name is in proper position: 18, 19 or 20)
        """
        old_name = resi.resname
        if len(old_name.strip()) == 1:
            if self.parameters.res_name == 20: return '  '+old_name.strip()
            elif self.parameters.res_name == 19: return ' '+old_name.strip() + ' '
            elif self.parameters.res_name ==18: return old_name.strip() + '  '
            else: return old_name
        else: return old_name
        
  
    def get_new_resi(self, resi):
        """
        Returns an empty Rediue object with proper name and flag.
        """
        new_id = self.get_new_id(resi)
        new_name = self.get_new_name(resi)
        new_resi = Residue(new_id, new_name, resi.segid)
        return new_resi
        
        
    def get_ribose_names_dict(self,  at_name_list):
        """
        Looks for all ribose atoms with wrong nomenclature in given residue.
        Returns dict {wrong_atom_name: name_that_should_be_there}
        """
        ribose_dict = {}
        # all ribose names should have *
        if self.parameters.ribose == '*' :
            for at in at_name_list:
                if "'" in at: ribose_dict[at] = at.replace("'", '*')
        # all ribose names should have '
        elif self.parameters.ribose =="'":
            for at in at_name_list:
                if "*" in at: ribose_dict[at] = at.replace("*", "'")
        return ribose_dict
        
        
    def get_p_names_dict(self, at_name_list):
        """
        Looks for all phosphate group atoms with wrong nomenclature in given residue.
        Returns dict {wrong_atom_name: name_that_should_be_there}
        """
        phosphate_dict = {}
        # all phosphate group atoms should have new nomenclature
        if self.parameters.phosphate == 'OP1':
            for at in ['O1P', 'O2P', 'O3P']:
                phosphate_dict[at] = 'OP' + at[1]                
        # all phosphate group atoms should have old nomenclature
        if self.parameters.phosphate == 'O1P':
            for at in ['OP1', 'OP2', 'OP3']:
                phosphate_dict[at] = 'O' + at[2] + 'P'                        
        return phosphate_dict
        
        
    def get_new_at_list(self, resi):
        """ 
        Returns complete list of atoms for residue.
        All atoms have proper nomenclature in respect of ribose and phosphate group.
        """
        at_list = [x for x in resi.child_list]
        if self.parameters.ribose == None and self.parameters.phosphate == None: return at_list
        at_name_list = [x.id for x in at_list]
        
        ribose_names = self.get_ribose_names_dict(at_name_list)
        p_names = self.get_p_names_dict(at_name_list)
        ribose_old_names = list(ribose_names.keys())
        p_old_names = list(p_names.keys())
        
        new_at_list = []
        for at in at_list:
            if at.id in ribose_old_names: new_at_list.append(self.change_at_name(at, ribose_names[at.id]))
            elif at.id in p_old_names: new_at_list.append(self.change_at_name(at, p_names[at.id]))
            # needs to be there because difference between new and old biopython
            else: new_at_list.append(self.change_at_name(at, at.id))
        return new_at_list
        
        
    def convert_resi(self, resi):
        """
        Converts residue.
        """
        # change resi flag and position of resi name
        new_resi = self.get_new_resi(resi)
        # change ribose and P group atom  names
        atom_list = self.get_new_at_list(resi)
        for at in atom_list:
            new_resi.add(at)
        return  new_resi
        
        
    def convert_resi_list(self,  resi_list):
        """
        Converts a list of residues.
        """
        new_resi_list = []
        for resi in resi_list:
            new_resi_list.append(self.convert_resi(resi))
        return new_resi_list
        
        
class ModernaStructureConverter(Converter):

    def __init__(self, st,  parameters):
        Converter.__init__(self, parameters)
        self.st = st # ModernaStructure obiect
        
    def convert_structure(self):
        """
        Converts ModernaStructure obiect.
        """
        new_resi = self.convert_resi_list([x for x in self.st])
        for resi in self.st:
            self.st.remove_residue(resi.identifier)
        for resi in new_resi:
            self.st.add_residue(resi)
            

    def save(self, output_file='converted.pdb'):
        """
        Writes ModernaStructure obiect with new nomenclature to a file.
        """
        self.st.write_pdb_file(output_file)
     
            
            
class PdbStructureConverter(Converter):
    
    def __init__(self, st,  parameters):
        Converter.__init__(self, parameters)
        self.st = st # PDB.Structure.Structure obiect       
        
        
    def convert_structure(self):
        """
        Converts PDB.Structure.Structure object.
        (all chains in all models)
        """
        for m in self.st:
            for ch in m:
                resi_list = []
                for resi in ch:
                    resi_list.append(resi)
                converted_resi = self.convert_resi_list(resi_list)
                for resi in resi_list:
                    ch.detach_child(resi.id)
                for resi in converted_resi:
                    ch.add(resi)
    

    def save(self,  output_file='converted.pdb'):
        """
        Saves structure to a file.
        """
        io = PDBIO()
        io.set_structure(self.st)
        io.save(output_file)        
        
        
def convert(input_file,  output_file='converted.pdb',  ribose=None, phosphate=None, res_name=None, flag=None):
    """Converts pdb file"""
    param = ConversionParameters(ribose, phosphate, res_name, flag)
    st = PDBParser().get_structure('', input_file)
    conv = PdbStructureConverter(st,  param)
    conv.convert_structure()
    conv.save(output_file)
