#!/usr/bin/env python
#
# ModernaAlphabet.py
#
# Defines alphabet for RNA including modifications.
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


from rna_tools.tools.mini_moderna3.moderna.Constants import MODIFICATION_NAMES_TABLE_PATH, ANY_RESIDUE, \
    RESIDUE_WITHOUT_ONE_LETTER_ABBREV 
from rna_tools.tools.mini_moderna3.moderna.util.Errors import AlphabetError


class AlphabetEntry(object):
    """
    Collects information about nomenclature for one residue.

    Arguments:
    - long abbreviation e.g. m1A, m66Am, A
    - short abbreviation - one letter abbreviation or ANY RESIDUE (X) 
        if there is no such abbreviation for residue
    - pdb abbreviation - abbreviation used in pdb e.g. MMA, MTA, UNK
    - full name
    - original base - stores name of original base for modification

    If a base is not modified then long_abbrev, short_abbrev, 
    pdb_abbrev and original_base are equal.
    """
    def __init__(self, long_abbrev=None, short_abbrev=None, pdb_abbrev=None, \
                 full_name=None, original_base=None, new_abbrev=None, \
                 category=None):
        self.long_abbrev = long_abbrev
        self.short_abbrev = short_abbrev
        self.pdb_abbrev = pdb_abbrev
        self.full_name = full_name
        self.original_base = original_base
        self.new_abbrev = new_abbrev
        self.category = category

    def __str__(self):
        return '<'+self.new_abbrev+'>'

    def __repr__(self):
        return self.new_abbrev
        
    def __eq__(self, other_ae):
        if self.short_abbrev == ANY_RESIDUE or other_ae.short_abbrev == ANY_RESIDUE: 
            return False
        elif self.long_abbrev == other_ae.long_abbrev: 
            return True
        else: 
            return False
        
    @property
    def purine(self):
        if self.original_base in ['A', 'G']: 
            return True
    
    @property
    def pyrimidine(self):
        if self.original_base in ['U', 'C']:
            return True

    def get_new_notation_tuple(self):
        """Converts an AlphabetEntry to a (string, length) tuple"""
        if self.short_abbrev in ['A', 'G', 'C', 'U', 'Q']:
            string = self.new_abbrev
            length = len(self.new_abbrev)
        else:
            string = '0' * (4 - len(self.new_abbrev)) + self.new_abbrev
            length = 4
        return string, length
        


class Alphabet(dict):
    """
Instance of this class contains a dict with nomenclature for modifications:
{ key: AlphabetEntry }
key - long abbreviation e.g. m1A, m66Am, A
AlphabetEntry conteins information as short abbreviation, pdb abbreviation, full name, original base.
To lern more see documentation for AlphabetEntry

there in an another dict available as Alphabet.one_letter_original:
{ key: AlphabetEntry}
key - short abbreviation (one letter code)
    """
    def __init__(self, table_path = MODIFICATION_NAMES_TABLE_PATH):
        dict.__init__(self) 
        self.table_path = table_path
        self.parse_alphabet_table()
        self._short_original = {}
        self._new_original = {}
        self.set_short_original()
        self.set_new_original()

    def parse_alphabet_table(self):
        """Parses table with different kinds of abbreviations."""  
        try:
            infile = open(self.table_path)
        except IOError:
            raise AlphabetError('File does not exist: %s ' % self.table_path) 

        for line in infile:
            if line.startswith('#') or line == '' or line == '\n':
                continue
            tokens = line.replace('\n','').split('\t')
            tokens = [t.strip() for t in tokens]
            if len(tokens) == 7: 
                aentry = AlphabetEntry()
                aentry.new_abbrev = tokens[0]
                aentry.original_base = tokens[1]
                aentry.long_abbrev = tokens[2]
                aentry.full_name = tokens[3]
                if len(tokens[4]) == 1:
                    aentry.short_abbrev = tokens[4]
                else:
                    aentry.short_abbrev = RESIDUE_WITHOUT_ONE_LETTER_ABBREV
                aentry.pdb_abbrev = tokens[5]
                aentry.category = tokens[6]
                self[tokens[2]] = aentry
            elif len(tokens) >0 :
                raise AlphabetError('Wrong line format: %s' %line)


    def set_short_original(self):
        """Creates short_original dict."""
        for abbrev in list(self.keys()):
            if self[abbrev].short_abbrev != RESIDUE_WITHOUT_ONE_LETTER_ABBREV:
                self._short_original[self[abbrev].short_abbrev] = self[abbrev]
        # add one defined entry for all entries without one letter abbreviations
        if RESIDUE_WITHOUT_ONE_LETTER_ABBREV in self:
            self._short_original[RESIDUE_WITHOUT_ONE_LETTER_ABBREV] = self[RESIDUE_WITHOUT_ONE_LETTER_ABBREV]


    def set_new_original(self):
        """Creates new_original dict."""
        for abbrev in list(self.keys()):         
            if self[abbrev].new_abbrev not in list(self._new_original.keys()):
                self._new_original[self[abbrev].new_abbrev] = self[abbrev]


    def get_short_original(self, short_abbrev):
        """Returns proper alphabet entry"""
        try:
            return self._short_original[short_abbrev]
        except KeyError:
            raise AlphabetError('This residue [%s] has no one letter abbreviation, cannot return alphabet entry'%short_abbrev)


    def get_new_original(self, new_abbrev):
        """Returns proper alphabet entry"""
        try:
            return self._new_original[new_abbrev]
        except KeyError:
            raise AlphabetError('This residue [%s] has no new abbreviation, cannot return alphabet entry'%new_abbrev)

# initialize the Alphabet once, for all other classes to import.
alphabet = Alphabet()
