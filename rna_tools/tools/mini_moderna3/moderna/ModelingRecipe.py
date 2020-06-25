#!/usr/bin/env python
#
# ModelingRecipe
#
# ModelingRecipe sorts alignment positions into different boxes (1 box == 1 operation, e.g. copying residues)
# 
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "rother.magdalena@gmail.com"
__status__ = "Production"

import os

from rna_tools.tools.mini_moderna3.moderna.sequence.RNAAlignment import RNAAlignmentParser, DEFAULT_SHRINK
from rna_tools.tools.mini_moderna3.moderna.Constants import STANDARD_BASES
from rna_tools.tools.mini_moderna3.moderna.util.Errors import AlignmentError
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log


MODES = ['has_gap', 'is_mismatch', 'has_template_gap', \
    'has_target_gap', 'is_unidentified']


class ModelingRecipe(object):
    """
    """
    def __init__(self):
        self.copy = []
        self.copy_backbone = []
        self.exchange = []
        self.remove_modifications = [] 
        self.add_modifications = []
        self.add_fragment = []
        self.difficult = [] # ??? report that an user should see this
        self.add_fragment_5p = []
        self.add_fragment_3p = []

        # MR: PROPOSAL: make orer of modeling_tasks flexible.
        #self.order = {'copy_backbone': 1, ...}
        # def set_order

class RecipeMaker(object):
    """
    All AlignmentPosition instances from an alignment are packed into boxes 
    which represent tasks for modeling.
    There are fallowing boxes:
    - copy ---> the same template and target letter and no ANY_RESIDUE
    - copy_backbone ---> target and template letter are present and one or two of them are ANY_RESIDUE
    - exchange ---> target and template letter are present and different and not ANY_RESIDUE
    - remove_modifications ---> template and target letters are present and template residue is not a standard one
    - add_modification ---> template and target letters are present and target residue is not a standard one
    - add_fragment ---> ther is no template or no target (never both) letter, this positions belongs to
    gap that doesn't have any other gap in target or template in distance less than 2
    - difficult ---> gaps different than above  
    """
    def __init__(self, alignment):
        self.recipe = ModelingRecipe()
        self.alignment = alignment
        self.set_alignment_properties()

    def get_differences (self, mode='is_different'):
        """
        Returns these positions in alignment which are different 
        in a given sequences as a list of AlignmentPosition instances.

        Arguments:
        - mode (by default 'all')
        Modes:
        - 'is_different' - all differences are taken into account 
            ('mismatch','gap in template','gap in target')
        - 'has_gap' - only gaps in template and gaps in target
        - 'has_template_gap'
        - 'has_target_gap'
        - 'is_mismatch'
        - 'is_unidentified' - when one of letters or both is ANY_RESIDUE
        """
        differences = []
        if mode in MODES:
            for apos in self.alignment:
                method = getattr(apos, mode)
                if method():
                    differences.append(apos)
            return differences
        raise AlignmentError("Bad mode type. Should be one of %s."%str(MODES))
        
    def remove_ap_from_box(self, apos):
        """
        Removes given AlignmentPosition from given boxes.
        By defult removes apos from all boxes besides add_fragment ons.
        """
        if apos in self.recipe.copy:
            self.recipe.copy.remove(apos)
        if apos in self.recipe.copy_backbone:
            self.recipe.copy_backbone.remove(apos)
        if apos in self.recipe.exchange:
            self.recipe.exchange.remove(apos)
        if apos in self.recipe.add_modifications:
            self.recipe.add_modifications.remove(apos)
        if apos in self.recipe.remove_modifications:
            self.recipe.remove_modifications.remove(apos)
        if apos in self.recipe.add_fragment:
            self.recipe.add_fragment.remove(apos)
        if apos in self.recipe.add_fragment_5p:
            self.recipe.add_fragment_5p.remove(apos)
        if apos in self.recipe.add_fragment_3p:
            self.recipe.add_fragment_3p.remove(apos)
        if apos in self.recipe.difficult:
            self.recipe.difficult.remove(apos)
        #TODO is this used?

    def segregate_gaps(self, gaps_devided):
        """ 
        Segregate given gaps list into proper boxes
        (add_fragmen, add_fragment_5p, add_fragment_3p, difficult).
        """
        for gap in gaps_devided: # gap is a tuple of ap belonging to one gap
            # ALL STRANGE CASES AT THE BEGINING AND AT THE END OF THE ALIGNMENT
            if gap[0].alignment_position == 1:
                # gap at beginning of alignment
                if gap[0].target_letter: 
                    self.recipe.add_fragment_5p.append(gap)

            elif gap[0].alignment_position == 2:
                # gap at the second position in the alignment
                ap_1before = self.alignment[gap[0].alignment_position -1]

                if gap[0].target_letter and ap_1before.target_letter: 
                # AA
                # A-
                    self.remove_ap_from_box(ap_1before)
                    self.recipe.add_fragment_5p.append([ap_1before]+ gap)

                elif gap[0].target_letter and ap_1before.template_letter:                       
                # -A
                # A-
                    self.recipe.add_fragment_5p.append(gap)

                elif ap_1before.target_letter and ap_1before.template_letter:
                # A-
                # AA
                    self.remove_ap_from_box(ap_1before)
                    self.recipe.add_fragment_5p.append([ap_1before])


            elif gap[-1].alignment_position == len(self.alignment):
                # gap at end of alignment
                if gap[0].target_letter: 
                    self.recipe.add_fragment_3p.append(gap)

            elif gap[-1].alignment_position == len(self.alignment)-1:           
                # gap at the second position before the end of the alignment
                ap_1after = self.alignment[len(self.alignment)]
    
                if gap[-1].target_letter and ap_1after.target_letter:
                # AA
                # -A
                    self.remove_ap_from_box(ap_1after)
                    self.recipe.add_fragment_3p.append(gap+[ap_1after])

                elif gap[-1].target_letter and ap_1after.template_letter:
                # A-
                # -A
                    self.recipe.add_fragment_3p.append(gap)

                elif ap_1after.target_letter and ap_1after.template_letter:
                # AA
                # -A
                    self.remove_ap_from_box(ap_1after)
                    self.recipe.add_fragment_3p.append([ap_1after])

            # REGULAR GAPS IN THE MIDDLE OF THE ALIGNMENT
            else:
                ap_1before = self.alignment[gap[0].alignment_position -1]
                ap_2before = self.alignment[gap[0].alignment_position -2]
                ap_1after = self.alignment[gap[-1].alignment_position +1]
                ap_2after = self.alignment[gap[-1].alignment_position +2]

                #if ap_1before.is_different('gap') or ap_2before.is_different('gap') \
                #or ap_1after.is_different('gap') or ap_2after.is_different('gap'):
                if ap_1before.has_gap() or ap_1after.has_gap():
                    self.recipe.difficult.append(gap)
                    #self.difficult.append(gap)
                    # MM CHANGES:
                elif  ap_2before.has_gap() or ap_2after.has_gap():
                    self.recipe.add_fragment.append([ap_1before] + gap + [ap_1after])
                    #TODO: too long gaps should remain difficult
                    #TODO: gaps close to each other should be inserted,
                    #      but with one base less cut on each side.
                    #      --> also 
                else: self.recipe.add_fragment.append([ap_2before, ap_1before] + gap + [ap_1after, ap_2after])
      

    def set_add_fragment_property(self, mode='has_template_gap'):
        """
        Fill in add_fragment an difficult propertis.
        """
        gaps = self.get_differences(mode) 
        if gaps:
            gaps2 = [ap.alignment_position for ap in gaps] # list of int
            gaps_devided = [] # nested list of int
            gaps2.sort()
            counter = gaps2[0]
            new_gap = []

            # group gaps
            for num in gaps2:
                if num == counter: # same gap
                    new_gap.append(self.alignment[num])
                    counter += 1
                else: # new gap
                    gaps_devided.append(new_gap)
                    new_gap = []
                    new_gap.append(self.alignment[num])
                    counter = num+1
            if new_gap: 
                gaps_devided.append(new_gap)            
            if gaps_devided: 
                self.segregate_gaps(gaps_devided)


    def set_alignment_properties(self):
        """
        Fill in alignment 'boxes'.
        See class description.
        """
        self.recipe.copy = self.alignment.get_identical_positions()
        self.recipe.copy_backbone = self.get_differences('is_unidentified')
        
        differences = self.get_differences('is_mismatch')
        self.recipe.exchange = [ap for ap in differences if ap.target_letter.long_abbrev in STANDARD_BASES and ap.template_letter.long_abbrev in STANDARD_BASES]
        self.recipe.add_modifications = [ap for ap in differences if ap.target_letter.long_abbrev not in STANDARD_BASES]
        self.recipe.remove_modifications = [ap for ap in differences if ap.target_letter.long_abbrev in STANDARD_BASES and ap.template_letter.long_abbrev not in STANDARD_BASES]

        self.set_add_fragment_property('has_template_gap')
        self.set_add_fragment_property('has_target_gap')
