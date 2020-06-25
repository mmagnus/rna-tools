#!/usr/bin/env python
#
# Isostericity.py
#  
# Module for modeling isosteric base pairs.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Pawel Skiba, Magdalena Musielak, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Pawel Skiba"
__email__ = "pw.skiba@gmail.com"
__status__ = "Prototype"

"""
Code Review: 2010-02-01 by KR.

The code is very well structured and clear to read.

It might be necessary to check the atom name issue in more detail.
(are always the atoms with same names superimposed on each other?)

I also added some small comments what we can look out for 
during some cleanup/refactoring phase, but I think
it is definitely too early for that.

Tests are 100% working here.

Keep up the good work!

look out for TODO's below.
"""


from rna_tools.tools.mini_moderna3.moderna.Constants import DATA_PATH
from IsostericityMatrices import IsostericityMatrices
from rna_tools.tools.mini_moderna3.moderna.RNAResidue import RNAResidue
from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.ModernaSuperimposer import ModernaSuperimposer
from rna_tools.tools.mini_moderna3.moderna.util.Errors import IsostericityError
from rna_tools.tools.mini_moderna3.moderna.analyze import BasePairCalculator
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log
import os

class Isostericity:
    """
    A class for modeling isosteric RNA base pairs
    """

    def __init__(self, old_bp=None, new_bp=None, cutoff=1.0, orientation="cis", force=False):
        """
        Arguments:
        - old_bp - single tuple/list: (moderna_residue1, moderna_residue2)
                   where moderna_residue is ModernaResidue object from model 
        - new_bp - new base pair description in format 'BP'
                   where BP is base pair shortcut
                   example: 'AC' - 'Adenine-Cytosine'
        - cutoff - maximum allowed RMSD
        - orientation - cis or trans
        - force  - change base pair even if there is no interaction between old_bp 
    
        Starts immediately if finds input data 
        """
        self.old_bp = old_bp
        self.new_bp = new_bp        
        self.cutoff = cutoff
        self.orientation = orientation[0].lower()
        self.force = force
        self.interaction = ''
        self.rmsd = None
        self.matrices = IsostericityMatrices()
        self.superimposer = ModernaSuperimposer()
        #self.pyrimidine_atoms = ["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "N1"]
       # self.purine_atoms = ["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "N9"]
        #self.pyrimidine_atoms = ["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "N*"]
        #self.purine_atoms = ["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "N*"]
        #self.pyrimidine_atoms = ["P", "O3'"]
        #self.purine_atoms = ["P", "O3'"]
        #self.pyrimidine_atoms = ["P", "O5'",  "C5'", "C4'", "C3'", "O3'"]
        #self.purine_atoms = ["P", "O5'",  "C5'", "C4'", "C3'", "O3'"]
        self.pyrimidine_atoms = ["C1'", "N*"]
        self.purine_atoms = ["C1'", "N*"]
        self.result_bp = None #contains result base pair after modelling
        if None not in (self.old_bp, self.new_bp): self.check_input_data()


    def check_input_data(self):
        """
        Checks input data propriety
        """
        log.write_message("--------------------Isosteric base pairs modelling--------------------\n")
        if type(self.old_bp) is not tuple and type(self.old_bp) is not list:
            raise IsostericityError("Input data must be a tuple or list")
        elif len(self.old_bp) is not 2:
            raise IsostericityError("Input data length must be equal 2")
        elif type(self.old_bp[0]) is not RNAResidue or type(self.old_bp[1]) is not RNAResidue:
            raise IsostericityError("First input element must be a tuple with two interacting ModernaResidue objects")
        elif type(self.new_bp) is not str or len(self.new_bp) is not 2:
            raise IsostericityError("Second input element must be a 2-chars String description of new base pair") 
        
        self.interaction = BasePairCalculator.base_pair_calc(self.old_bp[0], self.old_bp[1])
        
        if self.interaction is not None:
            self.interaction = self.orientation + self.interaction #BasePairCalculator doesn't return cis/trans            

        if self.interaction is None:
            if self.force == False:
                raise IsostericityError("There is no interaction between input residues")
            else:
                log.write_message("No interaction between input residues - force superposition")
                self.force_superimpose(self.old_bp, self.new_bp)
        elif not os.path.exists(DATA_PATH+'base_pairs/'+self.interaction+'_'+self.new_bp+'.pdb'):
            self.interaction = self.interaction[0]+self.interaction[2]+self.interaction[1]            
            if not os.path.exists(DATA_PATH+'base_pairs/'+self.interaction+'_'+self.new_bp+'.pdb'):
                raise IsostericityError("New base pair or interaction '"+ self.interaction + "_" + self.new_bp +"' is unknown")
        
        self.superimpose()
        log.write_message("----------------End of isosteric base pairs modelling-----------------\n")


    def superimpose(self):
        """
        Checks isostericity and runs BackboneBuilder if pairs are not isosteric
        or superimpose pairs if are isosteric 
        """        
        old_bp = self.old_bp[0].resname
        old_bp += self.old_bp[1].resname

        isosteric = self.matrices.check_isostericity(old_bp, self.new_bp, self.interaction, self.cutoff)

        if not isosteric:
            if self.force:
                log.write_message("Base pairs " + str(self.old_bp[0].number) + " " + str(self.old_bp[1].number) + " " + old_bp + "->" + self.new_bp + " are not isosteric - force superposition")
                self.force_superimpose(self.old_bp, self.new_bp)
            else:
                log.write_message("Base pairs " + str(self.old_bp[0].number) + " " + str(self.old_bp[1].number) + " " + old_bp + "->" + self.new_bp + " are not isosteric - starting backbone builder")
                self.run_backbone_builder()
        else:
            log.write_message("Base pairs " + str(self.old_bp[0].number) + " " + str(self.old_bp[1].number) + " " + old_bp + "->" + self.new_bp + " are isosteric - starting superposition")
            self.superimpose_bp()        
           

    def create_atoms_lists(self, moved_residues):
        """
        Returns moved, fixed and to_move atoms lists
        """                
        to_move = moved_residues['1'].child_list + moved_residues['2'].child_list
        fixed = self.old_bp[0].child_list + self.old_bp[1].child_list

        type_moved = (moved_residues['1'].get_resname(), moved_residues['2'].get_resname())
        type_fixed = (self.old_bp[0].get_resname(), self.old_bp[1].get_resname)

        moved = self.select_atoms(moved_residues['1']) + self.select_atoms(moved_residues['2'])
        fixed = self.select_atoms(self.old_bp[0]) + self.select_atoms(self.old_bp[1])

        return moved, fixed, to_move        
                 
    
    def select_atoms(self, residue):
        """
        Returns list with atoms for Superimposer
        """    
        if residue.resname in ('A','G'):
            atom_list = self.purine_atoms
        else:
            atom_list = self.pyrimidine_atoms

        atoms = []        
        for atom in atom_list:
            atoms.append(residue[atom])

        return atoms


    def superimpose_bp(self):
        """
        Superimpose two base pairs
        Uses self.moved_atoms and self.standard_atoms
        """
        filename = self.interaction+'_'+self.new_bp+'.pdb'
        log.write_message("File used for superposition: "+filename)

        moved_residues = ModernaStructure('file', DATA_PATH+'base_pairs/'+filename).moderna_residues
        moved, fixed, to_move = self.create_atoms_lists(moved_residues)

        self.rmsd = ModernaSuperimposer(fixed, moved, to_move).rmsd

        log.write_message("Moved atoms: "+str(moved))
        log.write_message("Fixed atoms: "+str(fixed))
        log.write_message("RMSD: "+str(self.rmsd))

        if self.rmsd <= self.cutoff:
            log.write_message("\nRMSD: " + str(self.rmsd) + "\n")
            self.assign_residue_numbers(moved_residues['1'], moved_residues['2'])
        else:
            log.write_message("\nRMSD bigger than offset " + str(self.rmsd) + ">" + str(self.cutoff) + " - starting backbone builder\n")
            self.run_backbone_builder()        
                

    def force_superimpose(self, old_bp, new_bp):
        """
        Superimposes base pairs if there is no interaction between old_bp and force mode is True
        """
        min_rmsd = 100.0
        pdb_name = ''
        base_pair = None
        path = DATA_PATH+'base_pairs/'

        for pdb in os.listdir(path):
            if new_bp in pdb:
                moved_residues = ModernaStructure('file', path+pdb).moderna_residues
                moved, fixed, to_move = self.create_atoms_lists(moved_residues)
                rmsd = ModernaSuperimposer(fixed, moved, to_move).rmsd
                if rmsd < min_rmsd:
                    pdb_name = pdb      
                    min_rmsd = rmsd
                    base_pair = moved_residues

        log.write_message("The best RMSD ("+str(min_rmsd)+") with structure: "+pdb_name)  
        if min_rmsd <= self.cutoff:
            log.write_message("RMSD <= cutoff (" + str(self.cutoff) + ")")
            self.assign_residue_numbers(base_pair['1'], base_pair['2'])
        else:
            log.write_message("RMSD > cutoff (" + str(self.cutoff) + ") - starting backbone builder")
            self.run_backbone_builder()


    def assign_residue_numbers(self, res1, res2):
        """
        Compares old and new residues data and change new if necessary
        """
        if res1.number is not self.old_bp[0].number:
            res1.change_number(str(self.old_bp[0].number))
        if res2.number is not self.old_bp[1].number:
            res2.change_number(str(self.old_bp[1].number))        
        self.result_bp = (res1, res2)


    def run_backbone_builder(self):
        """
        Runs backbone builder
        """
        pass


###############################################

if __name__ == '__main__':
    struct = ModernaStructure('file', '../test/test_data/rna_structures/1EHZ.pdb').moderna_residues

    test1 = Isostericity((struct['14'], struct['13']), 'AC').result_bp 
    test2 = Isostericity((struct['12'], struct['23']), 'CG').result_bp  
    test3 = Isostericity((struct['20'], struct['22']), 'CU').result_bp 
    test4 = Isostericity((struct['8'], struct['24']), 'AU', 1.5, True).result_bp 
