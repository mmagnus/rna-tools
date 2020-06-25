#!/usr/bin/env python
#
# CheckPdb.py
#
# Supports work with PDB.Structure objects.
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


import re,os

from Bio.PDB.Atom import Atom
from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.builder.PhosphateBuilder import TerminalPhosphateBuilder
from rna_tools.tools.mini_moderna3.moderna.analyze.ChainConnectivity import are_residues_connected, is_chain_continuous

from rna_tools.tools.mini_moderna3.moderna.Constants import MISSING_RESIDUE, UNKNOWN_RESIDUE_SHORT, PHOSPHATE_GROUP, RIBOSE, BACKBONE_ATOMS, AA_ATOMS, \
    BACKBONE_RIBOSE_ATOMS,  PHOSPHORYLATED_NUCLEOTIDES
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log


# what about aa

HEADER = """________________________________________________________________
        
    ModeRNA (C) 2009 by Magdalena Musielak & Kristian Rother, 
                         Tomasz Puton, Janusz M. Bujnicki 
    
    Check PDB file for any features that may prevent the file from being a ModeRNA template.
    
    support: lenam@amu.edu.pl
________________________________________________________________       

        \n"""

MESSAGE = """OVERALL INFORMATION:

The most common reason that causes problems while modeling RNA with ModeRNA 
is that a template is not compatible with the alignment. Possible reasons include:

1. Backbone discontinuity.

   Symptom: 
   There are underscore symbols in the sequence.
   e.g. structure.get_sequence() ---> 'A_AAAA'
   
   Explanation:
   It means that a structures sequence is not continous and can be caused by 
   some missing backbone atoms, nonproper bond lengths etc. It may also mean 
   that residues are simply not connected, e.g. when there are two blunt ends of a helix.
   
   You should include the _ in your alignment explicitly e.g.:
   A_AAAA
   A-AGAA
   ModeRNA can try to fix broken backbone with ModeRNA, which helps in many cases.
   See the fix_backbone() function for more details.


2. Single ligands, ions or water molecules.
   Symptom: 
   There are single underscore-dot symbols at the end of the sequence.
   e.g. structure.get_sequence() ---> 'AAAA_.'
   
   Explanation:
   It means there is one unidentified residue (ion) at the and of the structure, 
   that is not connected with the rest of the structure.
   
   This should be included in the alignment (see below) or the ions should be removed from
   the template structure using the clean_structure() function.
   AAAA_.
   AAAA--
   
   (see 3. if there are more such symbols).


3. Multiple water or ion molecules.

   Symptom: 
   There are multiple underscore-dot symbols at the end of the sequence.
   e.g. structure.get_sequence() ---> 'AAAA_.'
   e.g. structure.get_sequence() ---> 'AAAA_._._._._._._._._._.' (10 water molecules).
   
   Solution:
   The water chain should be removed with the clean_structure() function or included in the alignment:
   AAAA_._._._._._._._._._.
   AAAA--------------------  


4. Protein chain.
   
   Symptom:
   The entire chain has a sequence like ._._._. etc.
   
   Solution:
   Use a different chain. RNA can't be modeled from proteins yet.
   However, the clean_structure() function will remove them.


5. Missing or broken bases.
   Symptom: 
   There are single dots within the sequence.
   e.g. structure.get_sequence() ---> 'AA.AAAAA'
   
   Explanation:
   It means that one residue is not complete: a whole base or its part is missing
   and consequently the residue is marked in the alignment as an unknown residue.
   It is not possible to model another base in this place in the automatic mode. 
   Regardless of used alignment the residue will stay unidentified in the model
   (only backbone atoms will be copied to the model).
   To be on the safe side, the according residue should be removed manually.
   

All described above features need to be included in the alignment or the template structure need to be modified. 
One may check whether a PDB file (or alignment) needs to be cleaned up.
PDB file modification can be done with the clean_structure() function.


"""

class PdbController:
    """
    Finds structure irregularities and fixes them.
    """
    
    def __init__(self, ms):
        self.ms = ms
        self.water = [] 
        self.ions = []
        self.unidentified = []
        self.unidentified_rna = []
        self.rna_ligand = []
        self.phosphate = []
        self.P_trace = []
        self.P_missing = []
        self.aa = []
        self.disconnected_residues = []
        self.disordered = []
        self.standard = []
        #self.problematic_residues = [] # all residues from abov cathegories
        
        self.check_unidentified_residues()
        self.continuous = self.check_backbone_continuity()
        self.check_disconnected_residues()
        self.stars = self.check_stars()
        self.OP = self.check_OP()
        self.check_missing_P()


    def __str__(self):
        if not self.has_problems():
            return 'Chain OK'
        else:
            pdb_string='The chain contains some features that may cause problems during modeling:\n'
            if self.water: pdb_string+='- Water residues present in the chain: %s\n'%', '.join([x.identifier for x in self.water])
            if self.ions: pdb_string+='- Ion residues present in the chain: %s\n'%', '.join([x.identifier for x in self.ions])
            if self.aa: pdb_string+='- Aminoacids present in the chain: %s\n'%', '.join([x.identifier for x in self.aa])
            if self.unidentified_rna: pdb_string+='- Unidentified RNA residues present in the chain: %s\n'%', '.join([x.identifier for x in self.unidentified_rna])
            if self.P_trace: pdb_string+='- P trace present in the chain: %s\n'%', '.join([x.identifier for x in self.P_trace])
            if self.P_missing: pdb_string+='- Atom P missing in residue: %s\n'%', '.join([x.identifier for x in self.P_missing])
            if self.unidentified: pdb_string+='- Unidentified residues present in the chain: %s\n'%', '.join([x.identifier for x in self.unidentified])
            if self.rna_ligand: pdb_string +='- RNA ligand present in the chain: %s\n'%', '.join([x.identifier for x in self.rna_ligand])
            if self.phosphate: pdb_string +='- Phosphate present in the chain: %s\n'%', '.join([x.identifier for x in self.phosphate])
            if self.stars: pdb_string+='- Some residues have * in ribose atom names\n'
            if self.OP: pdb_string+='- Some residues have O1P or O2P instead of OP1 and OP2.\n'
            if not self.continuous: pdb_string+='- Chain is discontinuous.\n'
            if self.disconnected_residues: pdb_string+='- Chain contains disconnected residues: %s.\n'%', '.join([x.identifier for x in self.disconnected_residues])
            if self.disordered: pdb_string +='- Atoms of some residues may have alternative location (only the first coordinates set taken): %s\n'%', '.join([x for x in self.disordered])
        return pdb_string
            
    def __repr__(self):
        return str(self)
        
    def has_problems(self):
        """
        Checks whether PdbControler found any problems in the structure.
        """
        if self.water: return True
        if self.ions: return True
        if self.unidentified: return True
        if self.unidentified_rna: return True
        if self.rna_ligand: return True
        if self.phosphate: return True
        if self.P_trace: return True
        if self.P_missing: return True
        if self.aa: return True
        if self.stars: return True
        if self.OP: return True
        if not self.continuous: return True
        if self.disordered: return True
        return False


    def check_all_bb_missing(self,  resi):
        """
        Returns True when residue does not have all backbone and ribose atoms
        (BACKBONE_RIBOSE_ATOMS)
        """
        for at_name in BACKBONE_RIBOSE_ATOMS:
            if resi.has_id(at_name): return False
            #TODO: resi.is_backbone_complete ??
        return True


    def is_rna_ligand(self,  resi):
        """
        Returns True when residue is a RNA ligabd (e.g. ATP, base)
        """
        if resi.new_abbrev == 'ligand': return True
        if resi.alphabet_entry.category == 'standard RNA': return False
        if self.check_all_bb_missing(resi): return True
        if resi.long_abbrev in PHOSPHORYLATED_NUCLEOTIDES: return True
        return False


    def is_phosphate(self, resi):
        """
        Checks whether resi has less than 5 atoms and contain P, OP1 and OP2
        """
        if len(resi.child_list) > 5: return False
        if resi.has_id('P') and (resi.has_id('OP1') or resi.has_id('O1P')) and (resi.has_id('O2P') or resi.has_id('OP2')):
            return True
        return False
    
    
    def check_unidentified_residues(self):
        """
        Finds unidentified residues in the structure.
        Sorts them into water, ions, unidentified.
        """
        for resi in self.ms:
            if resi.long_abbrev == UNKNOWN_RESIDUE_SHORT:   
                if resi.child_list[0].name == 'O' and len(resi.child_list)==1: 
                    self.water.append(resi)
                elif resi.child_list[0].name == 'P' and len(resi.child_list)==1: 
                    self.P_trace.append(resi)
                else:
                    rna = True
                    for at in BACKBONE_ATOMS:
                        if not resi.has_id(at): rna = False
                    if rna: self.unidentified_rna.append(resi)
                    else:
                        aa = True
                        for at in AA_ATOMS:
                            if not resi.has_id(at): aa = False
                        if aa: self.aa.append(resi)
                        elif self.is_phosphate(resi): self.phosphate.append(resi)
                        elif len(resi.child_list) < 7: self.ions.append(resi)
                        else: self.unidentified.append(resi)
            elif self.is_rna_ligand(resi): self.rna_ligand.append(resi)
            else: 
                if resi.has_double_coord: self.disordered.append(resi.identifier)
                self.standard.append(resi)


    def check_stars(self):
        """
        Checks whether ribos atoms contain stars in the atom names.
        Returns True or False.
        """
        for resi in self.standard:
            if re.findall('\*', ''.join([at.name for at in resi.child_list])): return True
        return False


    def check_OP(self):
        """
        Checks whether structure contains irregular names foe oxygens bonded to phospate.
        Returns True or False.
        """
        for resi in self.standard:
            if self.is_typical_category_residue(resi):
                if re.findall('O1P', ''.join([at.name for at in resi.child_list])): return True
                if re.findall('O2P', ''.join([at.name for at in resi.child_list])): return True
                if re.findall('O3P', ''.join([at.name for at in resi.child_list])): return True
        return False
    
    
    def check_missing_P(self):
        """
        Checks whether any residue has no P atom.
        Returns True or False 
        """
        result = False
        for resi in self.standard:
            if self.is_typical_category_residue(resi):
                if not resi.has_id('P'): 
                    self.P_missing.append(resi)
                    result = True
        return result


    def check_backbone_continuity(self):
        """
        Checks whether chain is continous.
        Returns True or False'
        Takes only regular residues into account
        (water, ions and unidentified residues are excluded.) 
        """
        st = ModernaStructure('residues', self.standard + self.unidentified_rna)
        return is_chain_continuous(st)


    def is_typical_category_residue(self,  res):
        """
        Returns True when residue belongs to standard RNA, standard DNA or modified category
        False when unknown, ligand, synthetic, stereoisomer, insertion, missing, ' '
        """
        if res.category in ['standard RNA', 'standard DNA', 'modified']: return True
        else: return False


    def check_disconnected_residues(self):
        """
        Looks for residues that belong to the chain, are not connected to the chain 
        and represents synthetic residues, ligands or stereoisomers.
        """
        temp = []
        st = ModernaStructure('residues', self.standard+self.unidentified_rna)
        all_resi = [x for x in st]
        if len(all_resi) == 1: temp.append(all_resi[0].identifier)
        elif not self.continuous:
            for x in range(len(all_resi)-1):
                if not are_residues_connected(all_resi[x],  all_resi[x+1]):
                    if len(all_resi) == 2: temp += [y.identifier for y in all_resi]
                    elif x+1==len(all_resi)-1: temp.append(all_resi[x+1].identifier)
                    elif not are_residues_connected(all_resi[x+1],  all_resi[x+2]): temp.append(all_resi[x+1].identifier)

        # checking whether disconnected residues are not sandard or modified ones
        for resi in temp:
            if not self.is_typical_category_residue(st[resi]):
                self.disconnected_residues.append(st[resi])

#    def check_problematic_residues(self):
#        """
#        Goes through all lists (self.water, self.ions ...) and
#        makes one list with all resides reported as problematic.
#        It can be run only after other checking functions.
#        """
#        resi_lists = [self.water, self.ions, self.unidentified, self.unidentified_rna, self.rna_ligand, self.phosphate, \
#                      self.P_trace, self.P_missing, self.aa, self.disconnected_residues, self.disordered]
#        for l in resi_lists:
#            for resi in 

    def write_log(self, file_name='structure.log'):
        """
        Writes information about structure irregularities to the file.
        """
        #TODO: this could be returned as a string, so we could print it in 'moderna.py -t <filename>'
        f = open(file_name,'w')
        f.write(HEADER)
        f.write(str(self))


    def remove_property(self,  property_list,  mode='residues'):
        """
        Removes required property from the structure
        (e.g. water, ions, unidentified residues)
        Mode - defines type of data present in the property list,
        it can be either 'residues' or 'identifiers'.
        """
        ms_id = [x.identifier for x in self.ms]
        prop_id = []

        for resi in property_list: 
            if  mode=='residues': 
                if resi.identifier in ms_id:
                    self.ms.remove_residue(resi.identifier)
                    prop_id.append(resi.identifier)
            elif mode == 'identifiers':  
                if resi in ms_id: 
                    self.ms.remove_residue(resi)
                    prop_id.append(resi)
        
        st_id = list(set(x.identifier for x in self.standard)-set(prop_id))
        standard_new = []
        for x in self.standard:
            if x in st_id: standard_new.append(x)
        self.standard = standard_new
                
        return []       
    
    def add_P_atom(self):
        """
        Adds missing P atoms (just when O5' and C5' are already there).
        """
        ms_id = [x.identifier for x in self.ms]
        new_P_missing = []
       
        for resi in self.P_missing:
            if resi.identifier in ms_id:
                add_P = True
                for at in ["C4'","C5'", "O5'"]:
                    if not resi.has_id(at):
                        add_P = False
                        log.write_message("Could not add P atom to residue %s. %s atom missing." %(resi.identifier, at))
                if add_P:
                    try:
                        tp = TerminalPhosphateBuilder(resi, self.ms)
                        tp.build_phosphate_group()
                        log.write_message("Phosphate group constructed on residue %s." %(resi.identifier))
                    except: 
                        log.write_message("Could not add P atom to residue %s.")
                        new_P_missing.append(resi)
                else: 
                    new_P_missing.append(resi)
        self.P_missing = new_P_missing
                

    def change_at_name(self, resi, old_name, new_name):
        """
        Changes single atom name.
        """
        at = resi[old_name]
        element = new_name[0]
        if element not in "CHONSP":
            element = 'C'
        new_at = Atom(new_name, at.coord, at.bfactor, at.occupancy, at.altloc, ' '+new_name, at.serial_number,  element=element)
        resi.detach_child(old_name)
        resi.add(new_at)
        #print new_at.coord, new_at.bfactor, new_at.occupancy, new_at.altloc, new_name, new_at.serial_number


    def fix_at_names(self, stars=True, OP=True):
        """
        Changes:
        - * ---> ' 
        - O1P ---> OP1
        - O2P ---> OP2
        - O3P ---> OP3
        """
        for resi in self.ms:
            if OP:
                if resi.has_id('O1P'): self.change_at_name(resi, 'O1P', 'OP1')
                if resi.has_id('O2P'): self.change_at_name(resi, 'O2P', 'OP2')
                if resi.has_id('O3P'): self.change_at_name(resi, 'O3P', 'OP3')
                self.OP = False
            if stars:
                for at in list(resi):
                    if '*' in at.name: 
                        self.change_at_name(resi, at.name, ''+at.name.replace('*',"'"))
                self.stars = False


    def clean_structure(self, rm_water=True, rm_ions=True, rm_aa=True, exchange_stars=True, exchange_OP=True, rm_unidentified=True, add_P=True,  rm_P_trace=True, rm_ligand=True, rm_phosphate=True,  rm_disconnected_resi=True, rm_unidentified_rna=False):
        """
        Fixes features that may causes some problems during modeling.
        """
        if rm_water: self.water = self.remove_property(self.water) 
        if rm_ions: self.ions = self.remove_property(self.ions) 
        if rm_unidentified: self.unidentified = self.remove_property(self.unidentified) 
        if rm_unidentified_rna: self.unidentified_rna = self.remove_property(self.unidentified_rna)
        if rm_ligand: self.rna_ligand = self.remove_property(self.rna_ligand) 
        if rm_phosphate: self.phosphate = self.remove_property(self.phosphate) 
        if rm_P_trace: self.P_trace = self.remove_property(self.P_trace) 
        if rm_aa: self.aa = self.remove_property(self.aa)
        if rm_disconnected_resi: 
            self.disconnected_residues = self.remove_property(self.disconnected_residues)#,  'identifiers')
            self.continuous = self.check_backbone_continuity()
        if exchange_stars and exchange_OP: self.fix_at_names()
        elif exchange_OP: self.fix_at_names(False, True)
        elif exchange_stars: self.fix_at_names(True, False)
        if add_P: self.add_P_atom()
        return self.ms


    def get_structure(self):
        """
        returns structure.
        """
        return self.ms


    def write_structure(self, file_name='fixed_structure.pdb'):
        """
        Writes structure to a pdb file.
        """
        self.ms.write_pdb_file(file_name)

   
