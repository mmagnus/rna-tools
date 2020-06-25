#!/usr/bin/env python
#
# Renumerator.py
#
# Generates identifiers for residues from ModernaFragment objects.
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


from rna_tools.tools.mini_moderna3.moderna.util.Errors import RenumeratorError

class Renumerator:
    """
    Class for generating number for residues renumeration.
    Useful in fragments of all kinds
    """
    def __init__ (self, residues_list, identifiers_list):
        if len(residues_list) != len(identifiers_list): 
            raise RenumeratorError('Cannot generate numbers - residue list and number list have different size.')
        self.residues_list = residues_list
        #KR: if this is the only place residue_list is used, we maybe won't need it.
        self.identifiers_list = identifiers_list
        
    def letter_generator(self,  first_num): 
        CHARACTERS = 'ABCDEFGHIJKLMNOPQRSTUWXYZ'
        for x in CHARACTERS:
            yield first_num + x

    def number_generator(self,  first_num):
        try: first_num = int(first_num)
        except: RenumeratorError('Number can not be generated, first number cannot be converted into int')
        while 1:
            first_num += 1
            yield str(first_num)

  
    def divide_identifiers(self):
        """
        divides given list of input identifiers
        e.g.
        [None, None, None, '2', '3'] ---> [[None, None, None], ['2', '3']]
        ['5', '6', None, None, None, '7', '8'] ---> [['5', '6'], [None, None, None], ['7', '8']]
        ['1', '2', None, None, None, '3', '4', '5', None, None, None, '6', '7'] ---> [['1', '2'], [None, None, None], ['3', '4', '5'], [None. None, None], ['6', '7']]
        """
        divided_list = []
        num = self.identifiers_list[0]
        temp_list = []
        for x in self.identifiers_list:
            if (x and num) or (not x and not num):
                temp_list.append(x)
            else:
                divided_list.append(temp_list)
                temp_list = [x]
            num = x
        if temp_list: divided_list.append(temp_list)
        return divided_list    
    
    
    def find_first_number(self,  before,  middle,  after):
        """
        KR what when the first number has already a letter?
        Should it generate an exception?
        """
        #KR: yes. a NumberingError, maybe?
        if not before and not after:
            return '1'
        elif before:
            return before[-1]
        elif after and before and len(middle)>=26: 
            return str(int(after[0])-int(before[-1])-1)
        elif after and len(middle)>=26: 
            first = int(after[0])
            if first > len(middle):
                return str(int(after[0])-len(middle)-1)
            else:
                raise RenumeratorError('Cannot generate numbers - not enough room for %i numbers before %s.'%(len(middle), after[0])) 
        else: 
            return str(int(after[0])-1)
    
    
    def get_generator(self, first_id, before, middle, after):
        """
        Returns identifiers generator.
        Checks whether it should be letters or numbers generator.
        In case of numbers checks whether the new numbering is possible
        """
        if len(middle) < 26: 
            return self.letter_generator(first_id)
        if before and after:
            if int(after[0])-int(before[-1]) > len(middle): 
                return self.number_generator(first_id)
            else:
                raise RenumeratorError('Cannot generate numbers - not enough room for %i numbers between %s and %s.'%(len(middle), before[-1], after[0])) 
        elif before:
            return self.number_generator(first_id)
        else:
            return self.number_generator(first_id)
    
    
    def prepare_identifiers(self,  before,  middle,  after):
        """
        middle - query list with None elements
        before - list with identifiers before query or None
        after - list with identifiers after query or None
        """
        fn = self.find_first_number(before, middle, after)
        id_generator = self.get_generator(fn, before, middle, after)
        identifiers = []
        for x in middle: identifiers.append(next(id_generator))
        return identifiers
    
    
    def get_identifiers_list(self):
        """
        Returns complete new identifiers list for given query identifiers list.
        e.g.
        [None, None, None] ---> ['1A', '1B', '1C']
        [None, None, None, '2', '3'] ---> ['1A', '1B', '1C', '2', '3']
        ['5', '6', None, None, None, '7', '8'] ---> ['5', '6', '6A', '6B', '6C', '7', '8']
        ['1', '2', None, None, None, '3', '4', '5', None, None, None, '6', '7'] ---> ['1', '2', '2A', '2B', '2C', '3', '4', '5', '5A', '5B', '5C', '6', '7']
        """
        identifiers_list = []
        divided_list = self.divide_identifiers()
        for x in range(len(divided_list)):
            if divided_list[x][0] == None:
                if x>0: 
                    before = divided_list[x-1]
                else: 
                    before = None
                if x<len(divided_list)-1: 
                    after = divided_list[x+1]
                else: 
                    after = None
                identifiers_list += self.prepare_identifiers(before, divided_list[x], after)
            else: identifiers_list += divided_list[x]
        return identifiers_list


def renumber_section(struc, start, end, new_start):
    """
    Renumbers a part of a structure 
    from start..end, new numbers starting from new_start.
    """
    length = end-start+1
    for old, new in zip(list(range(start, start + length)), list(range(new_start, new_start+length))):
        struc.renumber_residue(str(old), str(new))
