#!/usr/bin/env python
#
# ModernaStructure.py
#
# RNAStructure with some functions for editing.
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


import re
from rna_tools.tools.mini_moderna3.moderna.RNAResidue import RNAResidue
from rna_tools.tools.mini_moderna3.moderna.RNAChain import RNAChain
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence
from rna_tools.tools.mini_moderna3.moderna.analyze.BasePairCalculator import base_pair_calc
from rna_tools.tools.mini_moderna3.moderna.builder.BackboneBuilder import BackboneBuilder
from rna_tools.tools.mini_moderna3.moderna.modifications import modify_residue
from rna_tools.tools.mini_moderna3.moderna.util.Errors import ModernaStructureError
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log

#TODO: PDB structures 3ftm and 3fic are messed up. 
# Should ModeRNA do something reasonable with them anyway?
# They have residues with negative numbers.

class ModernaStructure(RNAChain):
    """RNAChain extended by some editing functionality"""

    def find_residues_in_range(self, start_id=None, stop_id=None):
        """Returns identifiers of residues in the given range"""
        #KR: doesnt work because needs to be permissive about non-existing residues.
        residues_in_range = []
        if not start_id and not stop_id:
            raise ModernaStructureError('Could not find residues in range. No range given')
        elif start_id:
            in_range = False
        elif stop_id and not start_id: 
            in_range = True

        for resi in self:            
            if resi.identifier == str(stop_id): 
                in_range = False
            if in_range:
                residues_in_range.append(resi.identifier) 
            if resi.identifier == str(start_id): 
                in_range = True
        return residues_in_range
    
    def find_residues_not_in_range(self, start_identifier=None, stop_identifier=None):
        """Returns list of residues except these in given range"""
        residues = []
        residues_in_range = self.find_residues_in_range(start_identifier, stop_identifier)
        
        for res in self:
            if res.identifier not in residues_in_range: 
                residues.append(res)
        return residues
    
    def add_residue(self, resi, number=None, strict=True):
        """
        Adds a new residue.

        Arguments:
        * residue (as Bio.PDB.Residue or ModernaResidue instance)
        * residue number (by default original residue number)
        """
        #TODO: should be in RNAChain.
        residue_number = number or resi.identifier
        if residue_number in self.moderna_residues: 
            if strict:
                raise ModernaStructureError('There is a residue with the number "%s" already'%residue_number)
            else:
                self.remove_residue(residue_number)
        
        #TODO: is cloning really necessary?
        resi = RNAResidue(resi, resi.alphabet_entry)
        resi.change_number(residue_number)
        self.moderna_residues[residue_number] = resi



    def renumber_residue(self, old_number, new_number):
        """
        Changes residue identifier.

        Arguments:
        * an old residue identifier
        * a new residue identifier
        """
        if type(new_number) != str or new_number[0] not in [str(x) for x in range(0, 10)] or len(new_number)>5 or len(re.findall('[a-zA-Z]', new_number))>1:
            raise ModernaStructureError('Cannot change residue %s number for %s. The new number need to be string less than 5 characters long with digit in first position and including at most 1 letter.' % (str(old_number),  str(new_number)))
        resi = self[old_number]
        self.remove_residue(old_number)
        self.add_residue(resi, new_number)    
        log.write_message('Residue %s ---> renumbered to %s' % (str(old_number),  str(new_number)))
    
    def renumber_chain(self, start_residue='1'):
        """
        Changes chain numeration. Provides continous chain numeration 
        starting from given number (by default 1)

        Arguments:
        * start number (by default 1)
        """
        try:
            num = int(start_residue)
        except ValueError: 
            raise ModernaStructureError('Cannot start numeration with %s, requires number' % str(start_residue)) 
        
        self.sort_residues()
        counter = num
        temp_resi_list = []
        for resi in self:
            temp_resi = resi
            self.remove_residue(resi.identifier)
            temp_resi.change_number(str(counter))
            temp_resi.id = (temp_resi.id[0], temp_resi.id[1], ' ')
            temp_resi_list.append(temp_resi)
            counter += 1

        for resi in temp_resi_list:
            self.add_residue(resi)
        self.sort_residues()
        
        log.write_message('Chain was renumbered. The first residue is now %s.' % str(start_residue))
        
    def change_sequence(self, new_seq,  start_id=None,  stop_id=None):
        """
        Changes whole sequence (or indicated fragment) into a given one.

        Arguments:
        * new sequence
        * start residue identifier
        * stop residue identifier
        (by default all sequence is changed)
        """
        if type(new_seq) != Sequence:
            new_seq = Sequence(new_seq)
        if not start_id:
            start_id = self.first_resi.identifier
        if not stop_id:
            stop_id = self.last_resi.identifier
        if len(self[start_id:stop_id]) != len(new_seq):
            raise ModernaStructureError('Bad sequence length. Sequence should have %s characters.' %len(self[start_id:stop_id]))
        for res, seq_letter in zip(self[start_id:stop_id], new_seq):
            modify_residue(res, seq_letter.long_abbrev)

    def check_letters_in_residue_numeration(self):
        """Checks whether identifiers of residues contain any letters."""
        for resi in self:
            if re.findall('[a-z,A-Z]', resi.identifier):
                log.write_message('Structure contains residues with letters in numeration eg.: %s.' %resi.identifier)
                return True
        return False

    def get_modified_residues(self):
        """
        Returns a dict with all modified residues in the structure.

        key - residue identifier
        value - ModernaResidue instance
        """
        log.write_message('Looking for modifications.')
        modified_residues = {}
        for resi in self:
            if resi.modified:
                modified_residues[resi.identifier] = resi
                log.write_message('Residue %s: is modified (%s)'%(resi.identifier, resi.long_abbrev))
        if not modified_residues:
            log.write_message('There are no modifications in the structure.')
        return modified_residues

    def get_base_pairs(self, mode='dictionary'):
        """
        Returns a dictionary with all base pairs interactions in moderna structure - mode='dictionary'
        or yields a generator containing single residue with all interacting residues - mode='generator'
        """      
        #TODO: refactor
        def generator():
            was = []
            #TODO: use KDTree
            for resi1 in list(self):
                for resi2 in list(self):
                    if resi1 != resi2:
                        interaction = base_pair_calc(resi1, resi2)
                        if interaction != None and (resi1.identifier, resi2.identifier) not in was:
                            was.append((resi2.identifier, resi1.identifier))
                            yield interaction
        
        if mode == 'dictionary':
            interactions = {}
            for resi1 in list(self):
                interactions[resi1.identifier] = [base_pair_calc(resi1, resi2) \
                        for resi2 in list(self) if resi1 != resi2 and base_pair_calc(resi1, resi2) != None]
                        #TODO: now base_pair_calc is called twice!
            return interactions
        elif mode == 'generator':
            return generator()  
        else:
            return "Unknown mode selected: %s" %mode
            
    def contains_pseudoknot(self, a1, a2, bplist):
        """Checks if the pair (a1,a2) makes a pseudoknot with the list of base pairs. Returns boolean"""
        for b1, b2 in bplist:
            if (b1 < a1 < b2 < a2) or (b2 < a1 < b1 < a2) or (b2 < a2 < b1 < a1):
                return True
        return False
            
    def get_secstruc(self):
        """Returns secondary structure as a dot-bracket string."""
        #TODO: refactor out.
        bp_dict = self.get_base_pairs()
        res_indices = dict([(res.identifier, i) for i, res in enumerate(self)])
        ss = ['.']*len(self)
        bplist = []
        for resi in self:
            if resi.identifier not in bp_dict: continue
            res1 = resi.identifier
            for base_pair in bp_dict[res1]:
                # set bracket symbols.
                if base_pair.canonical or base_pair.wobble:
                    i1 = res_indices[res1]
                    i2 = res_indices[base_pair.resi2.identifier]
                    if i1 > i2:
                        i1, i2 = i2, i1
                    if not self.contains_pseudoknot(i1, i2, bplist):
                        ss[i1] = '('
                        ss[i2] = ')'
                        bplist.append((i1, i2))
        return ''.join(ss)
        
        
    def write_secstruc(self,  file_name='secstruc.vienna'):
        """Writes secstruc to a vienna file"""
        outfile = open(file_name,  'w')
        ss_str = '> Secondary structure:\n%s\n%s' % (self.get_sequence().seq_with_modifications.replace('_',''),  self.get_secstruc())
        outfile.write(ss_str)
        
    # ------------- helper functions --------------------
    def get_resi_before(self, resi):
        """Returns residue before the given residue"""
        resis = list(self)
        index = resis.index(resi)
        if index > 0:
            return resis[index-1]
            
    def get_resi_after(self, resi):
        """Returns residue after the given residue"""
        resis = list(self)
        index = resis.index(resi)
        if index < len(resis)-1:
            return resis[index+1]

    # --------------- BACKBONE REPAIR FUNCTION -------------------------
    def fix_backbone_between_resis(self, resi1, resi2):
        """Attempts to fix the backbone between the two given residues."""
        BackboneBuilder(resi1, resi2, self)
        
    def fix_backbone_after_resi(self, resi):
        """Fixes backbone on 3' side of a residue"""
        after = self.get_resi_after(resi)
        if after != None:
            self.fix_backbone_between_resis(resi, after)
        
    def fix_backbone_before_resi(self, resi):
        """Fixes backbone on 5' sides of a residue"""
        before = self.get_resi_before(resi)
        if before != None:
            self.fix_backbone_between_resis(before, resi)

    def fix_backbone(self):
        """Attempts to fix all backbone breaks."""
        log.write_message('\nChecking structure for interrupted backbones:')
        for i in range(2): # run twice to improve result
            #TODO: create better method for backbone fixing.
            last = None
            for resi in self:
                if last:
                    self.fix_backbone_between_resis(last, resi)
                last = resi
        log.write_message('Checking and repairing backbones finished.')
