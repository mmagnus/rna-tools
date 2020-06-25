#!/usr/bin/env python
#
# AlignmentMatcher.py
#
# Fixes Alignment object so that it matches a given template.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@genesilico.pl"
__status__ = "Production"


from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaAlphabet import alphabet
from rna_tools.tools.mini_moderna3.moderna.Constants import ANY_RESIDUE
from rna_tools.tools.mini_moderna3.moderna.util.Errors import AlignmentError

BREAK = alphabet['_']
GAP = alphabet['-']


class PairQueue(object):
    """
    Double queue algorithm that produces pairs of 
    AlignmentPosition and AlphabetEntry
    """
    def __init__(self, align, seq):
        # queues
        self.ap_queue = list(align)
        self.guide_queue = list(seq)
        self.ap_queue.reverse()
        self.guide_queue.reverse()
        # counters
        self.i_guide = 0
        self.i_ali = 0
        
    def has_more(self):
        """True as long as both queues have elements."""
        if self.ap_queue:
            return True
        
    @property    
    def pair(self):
        """Return current pair of elements."""
        guide = self.guide_queue[-1] if len(self.guide_queue) > 0 else GAP
        return guide, self.ap_queue[-1]
    
    def next_ap(self):
        """Moves one queue forward."""
        self.ap_queue.pop()
        self.i_ali += 1
        
    def next_guide(self):
        """Moves the other queue forward."""
        self.guide_queue.pop()
        self.i_guide += 1

    def next_both(self):
        """Moves both queues forward."""        
        self.next_guide()
        self.next_ap()

    
class AlignmentMatcher(object):
    """
    Compares a sequence to a sequence in an alignment,
    and fixes small differences by editing the alignment.
    """
    def __init__(self, alignment):
        self.align = alignment
        
    def is_template_identical(self, seq):
        """Returns boolean"""
        seq_t = seq.seq_without_breaks 
        seq_a = self.align.template_seq.seq_without_breaks
        return seq_t == seq_a
        
    def is_target_identical(self, seq):
        """Returns boolean"""
        return seq == self.align.target_seq
        
    def is_seq_fixable(self, seq):
        """Returns True if the guide sequence can be used for fixing."""
        if self.is_template_identical(seq): 
            return False
        tseq = self.align.template_seq
        if len(seq.seq_without_breaks) != len(tseq.seq_without_breaks):
            log.write_message("\nTemplate and alignment sequences differ in length - please check them manually.\n")
            return False
        return True

    def set_aligned_sequences(self, char_tuples):
        """Resets the sequences in the RNAAlignment object."""
        transposed = map(list, zip(*char_tuples))
        target = Sequence(transposed[0])
        template = Sequence(transposed[1])
        if len(target) != len(template):
            raise AlignmentError("Error correcting alignment; lenghts differ:\n%s\%s"%(str(target), str(template)))
        self.align.set_aligned_sequences(target, template)
        
    def check_breaks(self, guide, apos, dqueue, result):
        """Reacts on underscores in either of the sequences."""
        temp, targ = apos.template_letter, apos.target_letter
        if guide == BREAK and temp == BREAK:
            result.append((targ, temp))
            dqueue.next_both()
        elif guide == BREAK:
            log.write_message(".. break in template in position %i added to alignment."%(dqueue.i_guide+1))
            result.append((GAP, guide))
            dqueue.next_guide()
        else:
            log.write_message(".. break in alignment in position %i is not in template - ignored."%(dqueue.i_ali+1))
            dqueue.next_ap()

    def check_gaps(self, guide, apos, dqueue, result):
        """Reacts on gaps in the alignment."""
        if apos.has_template_gap():
            result.append((apos.target_letter, GAP))
        elif apos.has_target_gap():
            result.append((GAP, guide))
            dqueue.next_guide()
        dqueue.next_ap()
    
    def check_matches(self, guide, apos, dqueue, result):
        """Reacts on matches and mismatches."""
        temp, targ = apos.template_letter, apos.target_letter
        if temp.short_abbrev == ANY_RESIDUE:
            log.write_message(".. incomplete template residue in alignment position %i (%s/%s) - alignment edited." \
                              % (dqueue.i_ali+1, guide, temp))
            result.append((targ, guide))
        elif guide.short_abbrev == ANY_RESIDUE:
            log.write_message(".. unknown residue in alignment position %i (%s/%s) - alignment edited." \
                              % (dqueue.i_ali+1, guide, temp))
            result.append((targ, guide))
        elif guide.original_base != temp.original_base:
            log.write_message(".. different nucleobase in alignment position %i (%s/%s) - please check manually." \
                              % (dqueue.i_ali+1, guide, temp))
            result.append((targ, temp))
        elif guide != temp and guide.original_base == temp.original_base:
            log.write_message(".. different modified base found in alignment position %i (%s/%s) - alignment edited." \
                              % (dqueue.i_ali+1, guide, temp))
            result.append((targ, guide))
        elif guide == temp:
            result.append((targ, guide))
        else:
            # there may be cases not covered - report and ignore them
            log.write_message(".. don't know what to do about alignment position %i (%s/%s) - ignored." \
                              % (dqueue.i_ali+1, guide, temp))
            result.append((targ, temp))
        dqueue.next_both()

    def fix_template_seq(self, seq):
        """Adjusts the template sequence in the alignment to the given guide sequence."""
        # validate input seq
        if not self.is_seq_fixable(seq):
            return
        
        log.write_message("\nTemplate and alignment sequences differ - trying to fix small differences.\n")
        log.write_message("template           : %s"%seq)
        log.write_message("alignment (before) : %s\n"%\
                          self.align.aligned_template_seq)
        # iterate through positions
        dqueue = PairQueue(self.align, seq)
        result = []
        while dqueue.has_more():
            guide, apos = dqueue.pair
            if apos.has_gap():
                self.check_gaps(guide, apos, dqueue, result)
            elif guide == BREAK or apos.template_letter == BREAK:
                self.check_breaks(guide, apos, dqueue, result)
            else:
                self.check_matches(guide, apos, dqueue, result)

        self.set_aligned_sequences(result)
        log.write_message("\ntemplate          : %s"%seq)
        log.write_message("alignment (after) : %s\n"%str(self.align))

