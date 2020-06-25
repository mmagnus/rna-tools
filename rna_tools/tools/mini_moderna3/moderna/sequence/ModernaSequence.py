#!/usr/bin/env python
#
# ModernaSequence.py
#
# Enables work with RNA sequences.
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
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaAlphabet import alphabet
from rna_tools.tools.mini_moderna3.moderna.Constants import STANDARD_BASES
from rna_tools.tools.mini_moderna3.moderna.util.Errors import SequenceError

RE_NEW_NOMENCLATURE = re.compile('\d\d\d[A,C,U,G,a,c,t,g,X,N,H,<,;]')
SPECIAL_ABBREVS = ['a', 'c', 't', 'g', 'X', 'N', 'H', '<', ';', 'Q']

class Sequence(object):
    """
    Represents RNA sequences.
    
    Attributes:
    - seq_with_modifications ---> an original string given by a user
    - seq_without_modifications ---> a string in which all modified bases 
        are replaced by original ones 
    - has_modification 
        * False if all characters in given sequence are 
            A, C, G, U, ANY_RESIDUE or MISSING_RESIDUE
        * True if there are characters that represent modified residues.
    """

    def __init__(self, seq):
        """
        Arguments:
        -RNA sequence as a string
        """
        self._seq_without_modifications = ''
        self._seq_with_modifications = ''
        self._seq_new_notation = ''
        self.seq_alphabet_list = []
        self.has_modification = False      
        self.set_sequence(seq)

    def __str__(self):
        return self.seq_with_modifications

    def __repr__(self):
        return self.seq_with_modifications

    def __len__(self):
        return len(self.seq_alphabet_list)

    def __getitem__(self, args):
        """
        Allows to gain a character or characters from Sequence. 
        Returns AlpabetEntity object/list. 
        The starting point for counting is 0.
        """
        if type(args) == int:
            try: 
                return self.seq_alphabet_list[args] 
            except IndexError: 
                raise SequenceError('Sequence index out of range.')

        elif type(args) == slice:
            try: 
                return self.seq_alphabet_list[args.start:args.stop] 
            except IndexError: 
                raise SequenceError('Sequence index out of range.')

    def __iter__(self):
        return self.seq_alphabet_list.__iter__()

    def is_equal(self, seq1, seq2):
        """
        Checks if two RNA sequence objects are equal.
        If ANY_RESIDUE or MISSING_RESIDUE occurs in any of sequences 
        then in this position the sequence is equal to any other character. 
        """
        if len(seq1) != len(seq2):
            return False
        for char1, char2 in zip(seq1, seq2):
            if char1 != char2: 
                return False
        return True

    def __eq__(self, other):
        """
        Checks if two sequences are identical. Equal means that sequences length 
        is the same and the same bases on the same positions 
        however ANY_RESIDUE or MISSING_RESIDUE in one sequence is equal with any other character.

        Arguments:    
        - second sequence as a Sequence instance
        """
        if not self and not other:
            return True  
        elif self and other:
            if self.is_equal(self.seq_with_modifications, \
                        other.seq_with_modifications):
                return True
        return False

    def similar_to(self, other):
        """Checks if two sequences match, ignoring modifications."""
        if not self and not other:
            return True
        return self.seq_without_modifications == other.seq_without_modifications
 
    @property
    def seq_with_modifications(self):
        """Returns a string of the sequence with modified bases."""
        if not self._seq_with_modifications:
            self._seq_with_modifications = ''.join([e.short_abbrev for e in self])
        return self._seq_with_modifications
         
    @property
    def seq_without_modifications(self):
        """Prepares a sequence without modifications."""
        if not self._seq_without_modifications:
            self._seq_without_modifications = ''.join([e.original_base for e in self])
        return self._seq_without_modifications

    @property
    def seq_new_notation(self):
        """Returns a string with the sequence in JMB's new nomenclature."""
        if not self._seq_new_notation:
            result = ''
            for entry in self:
                if entry.new_abbrev in STANDARD_BASES:
                    result += entry.new_abbrev
                else:
                    result += '0' * (4-len(entry.new_abbrev)) + entry.new_abbrev 
            self._seq_new_notation = result
        return self._seq_new_notation
        
    @property
    def seq_without_breaks(self):
        """Returns sequence without break symbols."""
        return Sequence(str(self).replace('_', ''))

    @property
    def seq_without_gaps(self):
        """Returns sequence without gap symbols."""
        return Sequence(str(self).replace('-', ''))
 
    def set_sequence(self, seq):
        """
        Parses the sequence from a string 
        and produces a list of AlphabetEntries.
        """
        if type(seq) == str:
            x = 0
            while x < len(seq):
                if seq[x] == ' ': 
                    raise SequenceError('Sequence %s should not contain whitespace characters' % seq)
                elif seq[x] == '0':
                    if RE_NEW_NOMENCLATURE.match(seq[x:(x+4)]):
                        if seq[x+3] in SPECIAL_ABBREVS: 
                            new_abbrev = seq[x+3]
                        elif seq[x+1] == '0': 
                            new_abbrev = seq[x+2:(x+4)]
                        else:
                            new_abbrev = seq[x+1:(x+4)]                    
                        self.seq_alphabet_list.append(alphabet.get_new_original(new_abbrev))
                        x += 4
                    else:
                        raise SequenceError('Some irregularities in sequence: %s'  % seq[x:(x+4)])
                else:
                    self.seq_alphabet_list.append(alphabet.get_short_original(seq[x])) 
                    x += 1
        elif type(seq) == list: 
            self.seq_alphabet_list = seq
        else: 
            raise SequenceError('Bad argument type. Sequence instance takes string or list of alphabet_entities')
