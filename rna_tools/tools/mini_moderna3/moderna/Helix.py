#!/usr/bin/env python
#
# Helix
#
# Implements generation of A-type RNA helices.
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

from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.ModernaFragment import ModernaFragment53, ModernaFragment2D, \
    ModernaFragment553, ModernaFragment533
from rna_tools.tools.mini_moderna3.moderna.FragmentInsertion import FragmentInserter
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence
from rna_tools.tools.mini_moderna3.moderna.Renumerator import renumber_section
from rna_tools.tools.mini_moderna3.moderna.Constants import HELIX, HELIX_SUPERPOSITION, WC_BASE_PAIRS
from rna_tools.tools.mini_moderna3.moderna.util.Errors import ModernaFragmentError
from rna_tools.tools.mini_moderna3.moderna.analyze.ChainConnectivity import are_residues_connected


class Helix(ModernaStructure):
    """
    A-type RNA helix.
    """
    @property
    def strand5(self):
        """Return residues on the 5'strand."""
        return list(self)[:len(self)/2]

    @property
    def strand3(self):
        """Return residues on the 3'strand."""
        return list(self)[len(self)/2:]

    
class HelixBuilder(object):
    """
    Creates A-type RNA helices of a given length.
    """
    def __init__(self, data_type='file', data=HELIX, chain_name='A'):
        self.data_type = data_type
        self.data = data
        self.chain_name = chain_name
        
    def _get_helix(self):
        """Creates a helix structure."""
        helix = Helix(self.data_type, self.data, self.chain_name)
        if len(helix) % 2:
            raise ModernaFragmentError('Helix must have even length.')
        self.renumber_helix(helix)
        return helix
        
    def renumber_helix(self, helix):
        """
        Renumbers 5' strand starting from 1, and 3' strand starting from 401.
        """
        helix.renumber_chain()
        strand3 = helix.strand3
        first = int(strand3[0].identifier)
        last = int(strand3[-1].identifier)
        renumber_section(helix, first, last, 401)

    def make_helix_shorter(self, helix, dest_length):
        """Deletes base pairs from a helix."""
        resis = list(helix)[dest_length/2:-dest_length/2]
        for resi in resis:
            helix.remove_residue(resi.identifier)
        self.renumber_helix(helix)
        
    def make_helix_longer(self, helix): #, dest_length):
        """Adds base pairs to a helix."""
        helix2 = self._get_helix()
        helix2 = ModernaFragment53(helix2, anchor5=helix.strand5[-1], \
                                   anchor3=helix.strand3[0], strict=False)
        finsert = FragmentInserter()
        finsert.insert_fragment(helix2, helix)
        self.renumber_helix(helix)

    def build_helix(self, seq):
        """
        Returns a helix with the given sequence as a ModernaStructure object.
        """
        if not len(seq) % 2: 
            raise ModernaFragmentError('Expects odd number of characters \
in the sequence (FirstStrand_SecondStrand) e.g. "AAAA_GGGG"')
        helix = self._get_helix()
        dest_length = len(seq) - 1
        while len(helix) < dest_length:
            self.make_helix_longer(helix)
        if len(helix) > dest_length: 
            self.make_helix_shorter(helix, dest_length)
        if not len(helix) == dest_length:
            raise ModernaFragmentError('Helix construction failed. \
Length %i; should be %i'%(len(helix), dest_length))

        helix.change_sequence(seq.seq_without_breaks)
        return helix
    

class HelixFragmentBuilder(object):
    """
    Prepares helix fragments.
    Finds out which anchors on the upper section are present 
    and determines the according fragment type.
    """
    def _identify_upper_anchors(self, model, anchor5, anchor3):
        """Find residues connecting to anchors if any"""
        upper5, upper3 = None, None
        if model:
            upper = model.find_residues_in_range(anchor5.identifier, anchor3.identifier)
            if len(upper)>0:
                if are_residues_connected(anchor5, model[upper[0]]):
                    upper5 = model[upper[0]]
                if are_residues_connected(model[upper[-1]], anchor3):
                    upper3 = model[upper[-1]]
        return upper5, upper3

    def _add_anchors_to_seq(self, seq, upper5, upper3):
        """Eventually adds upper anchors to sequence."""
        # add lower anchors in any case
        seq = 'G'+str(seq)+'C'
        # determine upper anchors
        if upper5:
            upper5 = upper5.short_abbrev
            if not upper3:
                upper3 = WC_BASE_PAIRS[upper5]
        if upper3:
            upper3 = upper3.short_abbrev
            if not upper5:
                upper5 = WC_BASE_PAIRS[upper3]
        # add upper anchors
        if upper5 and upper3:
            strand_len = len(seq)//2
            seq = seq[:strand_len] \
                + upper5 + '_' + upper3 \
                + seq[-strand_len:]
        return Sequence(seq)
        
    def _create_frag(self, helix, seq, anchor5, anchor3, upper5, upper3, model):
        """Returns a ModernaFragment class"""
        seq = seq.seq_without_breaks
        if upper5 and upper3:
            return ModernaFragment2D(helix, \
                anchor5=anchor5, anchor3=anchor3, \
                anchor5_upper=upper5, \
                anchor3_upper=upper3, \
                frag5_upper=helix.strand5[-1], \
                frag3_upper=helix.strand3[0], \
                new_sequence=seq, \
                superposition_atoms=HELIX_SUPERPOSITION, 
                model=model
                )
        elif upper5:
            return ModernaFragment553
        elif upper3:
            return ModernaFragment533
        else:
            return ModernaFragment53(helix, \
                anchor5=anchor5, anchor3=anchor3, \
                new_sequence=seq, \
                sup_r5=HELIX_SUPERPOSITION, \
                sup_r3=HELIX_SUPERPOSITION, \
                strict=False
            )

    def build_fragment(self, anchor5=None, anchor3=None, \
                sequence=None, model=None):
        """Builds a Helix fragment."""
        upper5, upper3 = self._identify_upper_anchors(model, anchor5, anchor3)
        anchor_seq = self._add_anchors_to_seq(sequence, upper5, upper3)
        helix_builder = HelixBuilder()
        helix = helix_builder.build_helix(anchor_seq)
        return self._create_frag(helix, sequence, anchor5, anchor3, upper5, upper3, model)        

