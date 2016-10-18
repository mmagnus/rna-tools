#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The basic representation of RNA secondary structure.
"""

__author__ = "Tomasz Puton"
__credits__ = "Ewa Tkalińska, Łukasz Kozłowski, Kristian Rother"
__license__ = "GNU GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

import re, copy
from knotted2nested.rna2d import float_from_string
from itertools import cycle

class SecstrucError(Exception): pass
class UnknownElementError(Exception): pass
class PseudoknotTokensError(Exception): pass
class ConflictInBasePairsError(Exception): pass

class Secstruc:
    """
    Secondary structure string with indices for each position.
    Can handle pieces of secondary structure which are interrupted.
    
    Author: Kristian Rother
    """
    def __init__(self,secstr,indices=None):
        self.secstruc = secstr
        if indices == None:
            indices = range(len(self))
        if len(self.secstruc) != len(indices):
            raise SecstrucError("Cannot create Secstruc object (%s %s)"%(secstr,str(indices)))
        self.indices = indices

    def find_overhang5(self):
        """Returns overhang on the 5' end as secstruc objects."""
        n = len(self)
        for i in range(2,n):
            overhang = self[:i]
            if overhang.get_type() == "5'-overhang":
                return overhang

    def find_overhang3(self):
        """Returns overhang on the 3' end as secstruc objects."""
        n = len(self)
        for i in range(n-1,0,-1):
            overhang = self[i:]
            if overhang.get_type() == "3'-overhang":
                return overhang

    def find_helices(self):
        # return self.extract_elements('helix')
        helices = []
        begin = 0
        end = 0
        length = 0

        for i,j in self.get_base_pair_indices():
            # check if old helix continues
            if i == begin+length and j==end-length:
                length += 1
            else:
                # save old helix
                if length >0:
                    helices.append(self[begin:begin+length] + self[end-length+1:end+1])
                # new helix starts
                begin = i
                end = j
                length = 1
        if length > 0:
            helices.append(self[begin:begin+length] + self[end-length+1:end+1])
        return helices
 
    def get_base_pair_indices(self):
        """generates (i,j) tuples iterating through base pairs."""
        result = []
        left_bases = []
        n = len(self.secstruc)
        i = 0
        while i < n:
            if self.secstruc[i] == '(':
                left_bases.append(i)
            elif self.secstruc[i] == ')':
                left = left_bases.pop()
                result.append((left, i))
            i += 1
        result.sort()
        return result

    def get_bp_array(self):
        """
        Returns list:
        -1 : no pair
         n : index of paired base
        """
        a = [-1] * len(self.secstruc)
        for left, right in self.get_base_pair_indices():
            a[left] = right
            a[right] = left
        return a
 
    def find_loops(self):        
        return self.find_elements('loop')
    
    def find_nested_regions(self):
        """Helper function finding junctions and bulges 
        (things that have helices inside)."""
        for i,j in self.get_base_pair_indices():
            # check if i,j enclose a junction.
            nested_bps = 0
            skip_helix = 0
            junction = self[i]
            k = i+1
            while k < j:
                if self.secstruc[k] == ')':
                    skip_helix -= 1
                    if skip_helix == 0:
                        nested_bps += 1
                if skip_helix == 0:
                    junction += self[k]
                if self.secstruc[k] == '(':
                    skip_helix += 1
                k += 1
            junction += self[j]
            yield junction, nested_bps
 
    def find_junctions(self):
        """Returns a list of all junctions.
        Any region between two paired bases with more than 1 nested base pair
        Inside is a junction.
        """
        junctions = []
        for region, nested_bps in self.find_nested_regions():
            if nested_bps >= 2:
                junctions.append(region)
        return junctions
        
    def find_bulges(self):        
        """Returns a list of bulges."""
        bulges = []
        for region, nested_bps in self.find_nested_regions():            
            if nested_bps == 1 and str(region).find('.')>-1:
                bulges.append(region)
        return bulges


    def find_elements(self, element_type):
        """Finds elements of a given type."""
        result = []
        n = len(self)
        i = 0
        while i < n-1:
            j = n
            while j > i:
                ele = self[i:j]
                if ele.get_type() == element_type:
                    result.append(ele)
                j -= 1
            i += 1
        return result
                            

    def find_secstruc_elements(self):
        """
        Takes a secondary structure, and returns a list of
        component secondary structure elements.
        """
        result = []
        overhang = self.find_overhang5()
        if overhang: result.append(overhang)
        overhang = self.find_overhang3()
        if overhang: result.append(overhang)
        loops = self.find_loops()
        helices = self.find_helices()
        bulges = self.find_bulges()
        junctions = self.find_junctions()
        result += loops + bulges + helices + junctions
        return result

    def get_type(self):
        """Returns a type classification of this secondary structure as a string,
        if it is some basic type."""
        if re.search('^\(\(+\)\)+$',self.secstruc):
            # check both stems have equal length
            half = len(self.secstruc)/2
            if len(self.secstruc)%2 == 0 \
               and re.search('^\(+$',self.secstruc[:half]) \
               and re.search('^\)+$',self.secstruc[half:]):
                return 'helix'
            else: return None
        if re.search('^\.+\($',self.secstruc): return "5'-overhang"
        if re.search('^\)\.+$',self.secstruc): return "3'-overhang"
        if re.search('^\(\.\.+\)$',self.secstruc): return "loop"
        if re.search('^\(\(\)\.+\)|\(\.+\(\)\)$',self.secstruc): return "bulge"
        if re.search('^\((.*\(\)\.*)+\)$',self.secstruc): return 'junction'
        return None


    def __eq__(self, other):
        if not isinstance(other,Secstruc):
            other = Secstruc(other)
        if self.secstruc == other.secstruc \
           and self.indices == other.indices:
            return True
        
    def __add__(self, other):
        newsec = self.secstruc + other.secstruc
        newind = self.indices + other.indices
        return Secstruc(newsec, newind)

    def __getitem__(self,index):
        if type(index) == slice:
            newsec = self.secstruc[index.start:index.stop]
            newind = self.indices[index.start:index.stop]
            return Secstruc(newsec, newind)
        else:
            return Secstruc(self.secstruc[index],[self.indices[index]])

    def __repr__(self):
        return self.secstruc + ';' + str(self.indices)

    def __len__(self):
        return len(self.secstruc)

    def __str__(self):
        return self.secstruc

PSEUDOKNOTS_TOKENS = [
        ('[',']'), ('{','}'), ('<','>'), ('a','A'), ('b','B'), ('c','C'), \
        ('d','D'), ('e','E'), ('f','F'), 
        ('g','G'), ('h','H'), ('i','I'), ('j','J'), ('k','K'), ('l','L'), 
        ('m','M'), ('n','N'), ('o','O'), ('p','P'), ('q','Q'), ('r','R'), 
        ('s','S'), ('t','T'), ('u','U'), ('v','V'), ('w','W'), ('x','X'), 
        ('y','Y'), ('z','Z')]

## CLASSES Partners AND StructureString WERE COPIED FROM rna2d.py MODULE     ##
## FROM PyCogent WITHOUT ANY CHANGE.                                         ##
## (PyCogent ver. 1.3.1)                                                     ##

class Partners(list):
    """Holds list p such that p[i] is the tuple with indexes of the
        partners of i, or None.
    
        Primarily useful for testing whether a specified base is paired and,
        if so, extracting its partner(s).

        Each base may have 0 or more partners. If A pairs with B, B must pair
        with A.
    """
    
    def __setitem__(self, index, item):
        """Sets self[index] to item, enforcing integrity constraints."""
        if not type(index) == int:
            raise TypeError("index should be int, got %s" % str(type(index)))
            
        if item is not None and not isinstance(item, tuple):
            raise TypeError("item should be tuple, got %s" % str(type(item)))

        if item is None:
            
            previous = self[index]
            list.__setitem__(self, index, None)
            if previous is not None:
                for prev in previous:
                    # Trzeba sprawdzić powiązania i nie usuwać wszystkiego jak
                    # leci!
                    #import pdb; pdb.set_trace()
                    if len(self[prev]) == 1:
                        list.__setitem__(self, prev, None)
                    else:
                        decupled = tuple(sorted(list(set(self[prev]) - \
                                                     set([index]))))
                        list.__setitem__(self, prev, decupled)
                        
        else:
            if index in item:
                raise ValueError("Can't set base %s to pair with itself" % item)
            else:
                
                for i in item:
                    self_i_connected_to = None
                    if self[i] is None:
                        self_i_connected_to = self[index]
                        list.__setitem__(self, i, (index,)) # Tutaj tracimy
                                                            # info o usunięciu
                                                            # powiązania!
                        if self_i_connected_to is not None:
                            for sict in self_i_connected_to:
                                for sict_ in self[sict]:
                                    if sict not in self[sict_]:
                                        list.__setitem__(self, sict, None)
                        
                        if self[index] is None:
                            list.__setitem__(self, index, (i,))
                        else:
                            
                            curr_pos = self[index]
                            list.__setitem__(self, index, item)
                            if curr_pos == self[index]:
                                for cp in curr_pos:
                                    list.__setitem__(self, index,
                                                     tuple(sorted([i, cp])))
                    else:
                        list.__setitem__(self, i,
                                         tuple(sorted(list(self[i]) + [index])))
                
                   
                if index not in self[self[index][0]]:
                    list.__setitem__(self, index,
                                     tuple(sorted(list(self[i]) + [index])))
                    
                for i in item:
                    if i not in self[index]:
                        list.__setitem__(self, index,
                                         tuple(sorted(list(self[index]) + [i])))
                
                for index, partner in enumerate(self):
                    if partner is not None:
                        for p in partner:
                            if not index in self[p]:
                                list.__setitem__(self, index, None)
                                
                for index, partner in enumerate(self):
                    nonred = None
                    if partner is not None:
                        nonred = tuple(sorted(list(set(partner))))
                    list.__setitem__(self, index, nonred)
                
    def toPairs(self):
        """Converts the partners to sorted list of pairs."""
        result = set()
        for i, item in enumerate(self):
            if item is None:
                continue
            else:
                for i_ in item:
                    result.add(tuple(sorted([i, i_])))
                    
        return BasePairs(sorted(result))

    def _not_implemented(self, *args, **kwargs):
        """Raises NotImplementedError for 'naughty' methods.
        
        Not allowed any methods that insert/remove items or that change the
        order of the items, including things like sort or reverse.
        """
        raise NotImplementedError

class StructureString(str):
    """Base class for ViennaStructure and WussStructure. Immutable.
    
    StructureString holds a structure and a energy. By default energy is 
    set to None. 
    If you compare two StructureStrings the structure is the only important
    thing, since the energy is relative to the associated sequence.
    """
    Alphabet=None
    StartSymbols = {}      #dict of symbols that start base pairs
    EndSymbols = {}        #dict of symbols that end base pairs
    
    def __new__(cls, Structure, Energy=None):
        """Returns new StructureString."""
        a = cls.Alphabet
        if a:
            for i in Structure:
                if i not in a:
                    raise ValueError,\
                    "Tried to include unknown symbol '%s'" % i
        
        return str.__new__(cls,Structure)

    def __init__(self, Structure='', Energy=None):
        """Initializes StructureString with Structure and optionally Energy."""
        self.Energy = Energy
        self.toPartners()

    def __str__(self):
        """Returns string representaion of structure and energy, if known.

        Energy = 0 is different from Energy = None. Latter case is not printed.
        """
        if not self.Energy == None:
            return self + ' (' + str(self.Energy) + ')'
        else:
            return str.__str__(self)

    def toPartners(self):
        """Makes list containing partner of each position.
        
        Constructs a list from 0 to the number of bases, where each position
        contains the index of its pair (or None if it is unpaired).

        Note that the numbering starts at 0 for the first position.

        The algorithm here relies on the fact that any completely nested
        base-paired structure (no pseudoknots!) can be formally represented
        as a tree. Consequently, when you hit a closed base pair, you know
        that it it must pair with the last base pair you opened.
        """
        num_bases = len(self)
        result = [None] * len(self)
        stack = []
        start = self.StartSymbols
        end = self.EndSymbols
        for i, symbol in enumerate(self):
            if symbol in start:
                stack.append(i)
            elif symbol in end:
                curr = stack.pop()
                if result[i] is not None and result[curr] is not None:
                    result[i].append(curr)
                    result[curr].append(i)
                else:
                    result[i] = [curr]
                    result[curr] = [i]
               
        #test whether there are any open pairs left unaccounted for        
        if stack:
            raise IndexError, \
            "Too many open pairs in structure:\n%s" % self
        
        result = map(lambda x: tuple(x) if x is not None else None, result)
        return Partners(result)

    def toPairs(self, offset=1):
        """Makes list of (upstream,downstream) partners.
        
        Note that the numbering starts at 1 for the first position
        (not from 0, as usually in Python).
        Key will always be smaller than value.

        Result is in arbitrary order.
        
        offset: default is 1; offset is used to change the starting number
        """
        result = {}
        stack = []
        start = self.StartSymbols
        end = self.EndSymbols
        for i, symbol in enumerate(self):
            i_offset = i + offset
            if symbol in start.keys():       #open a pair
                stack.append((i_offset, symbol))
            if len(stack) > 0 and symbol in end.keys():
                for item in reversed(stack):
                    if item[1] == end[symbol]:
                        stack.remove(item)
                        result[i_offset] = item[0]
                        break
                    
        #test whether there are any open pairs left unaccounted for        
        if stack:
            raise IndexError, \
            "Too many open pairs in structure:\n%s" % self
        return BasePairs([(key, result[key]) for key in result]).directed()

##                                                                           ##
##                                                                           ##

class BasePairs(list):
    """The class was based on cogent.struct.rna2d.Pairs
    (PyCogent ver. 1.3.1). It has functionality necessary for processing 
    pseudoknots. 

    BasePairs should always be initialized in the following way:
    pairs = BasePairs(bplist),
    where bplist is a list with sub-lists containing base pair indices 
    (numbering starts from 1, not from 0, as usually in Python)!
    
    Adapted by: Tomasz Puton
    Email: t.puton@amu.edu.pl
    """
    def __init__(self, args):
        """Initializes BasePairs.
        
        args: a list of lists - each sub-list contains RNA base pair indices
        """
        list.__init__(self, args)
        self._pseudoknots = None
        self._conflicting = None
        
    def __hash__(self):
        """Let's make BasePair objects hashable, so they can be used as keys in
        dictionaries. E.g. for self equal [[2, 10], [4, 6], [3, 9], [5, 8]] 
        hash will be 210463958.
        """
        str_int = ""
        for base_pair in self:
            str_int += str(base_pair[0]) + str(base_pair[1])
        return int(str_int)
    
    def __cmp__(self, other):
        """Compares a hash value of self.directed() with a hash value of 
        other.directed()
        """
        directed_self = self.directed()
        try:
            directed_other = other.directed()
        except AttributeError:
            return 1
        else:
            return cmp(directed_self, directed_other)
    
    def directed(self):
        """Returns copy of self where all pairs are (upstream, downstream).
        
        Omits any unpaired bases and any duplicates. 
        Result is sorted according to base numbers (indices).
        """
        def cmp_for_bplist(x, y):
            if len(x) == 0 or len(y) == 0:
                return 0
            return cmp(x[0], y[0])
            
        seen = {}
        for up, down in self:
            if (up is None) or (down is None):
                continue    #omit unpaired bases
            if up > down:
                up, down = down, up
            seen[(up, down)] = True
        result = seen.keys()
        result.sort(cmp=cmp_for_bplist)
        return BasePairs(result)
    
    def symmetric(self):
        """Retruns copy of self where  each up, down pair has a down, up pair.
        
        Result is in arbitrary order. Double pairs and pairs containing None 
        are left out.
        """
        result = self.directed()
        result.extend([(down, up) for up, down in result])
        return BasePairs(result)

    def identifyNestedBasepairs(self):
        """Finds nested base pairs in bplist e.g. [2,21] is nested with respect
        to [1, 22]. Returns dict, where key is a base pair representative for
        a set of nested base pairs, while dict value associated with it is a
        list of nested base pairs.
        
        IMPORTANT! Base pair is considered nested if it does not form a 
        psuedoknot with another base pair. Current implementation sometimes
        fails to do that, so use this method with caution.
        
        bplist: a list of lists - every sub-list contains two integers which 
        are base numbers (numbering starts from 1!!! NOT from 0 as usually 
        in Python).
        """
        nested = {}
        sorted_bplist = self.directed()
        for base_pair in sorted_bplist:
            # base_pair is a tuple (thanks to .directed() method)
            if len(nested.keys()) == 0:
                # this is why we encapsulate base_pair in list and then feed
                # it to BasePairs
                nested[BasePairs([base_pair])] = []
            else:
                # we will always append base_pairs to the heighest key from
                # nested.keys() - this is why we have to sort it frist!
                last_key = sorted(nested.keys())[-1]
                if last_key[0][0] < base_pair[0] \
                and last_key[0][1] > base_pair[1]:
                    list_of_nested_base_pairs_of_last_key = \
                        nested[sorted(nested.keys())[-1]]
                    if len(list_of_nested_base_pairs_of_last_key) == 0:
                        nested[last_key].append(base_pair)
                    else:
                        # let's check whether base pair being analysed 
                        # is nested with respect to the last base_pair
                        if list_of_nested_base_pairs_of_last_key[-1][0] \
                        < base_pair[0] \
                        and list_of_nested_base_pairs_of_last_key[-1][1] \
                        > base_pair[1]:
                            nested[last_key].append(base_pair)
                        else:
                            nested[BasePairs([base_pair])] = []
                elif last_key[0][0] < base_pair[0] \
                and last_key[0][1] < base_pair[1]:
                    nested[BasePairs([base_pair])] = []            
        return nested

    def hasConflicts(self):
        """Returns True if the pair list contains conflicts.

        Conflict occurs when a base has two different partners, or is asserted
        to be both paired and unpaired.
        """
        self._conflicting = tuple()
        partners = {}
        for first, second in self:
            if first is None:
                if second is None:
                    continue    #no pairing info
                else:
                    first, second = second, first   #swap order so None is 2nd
            if second is None: #check first isn't paired
                if partners.get(first, None) is not None:
                    self._conflicting = tuple(sorted((first, second)))
                    return True
                else:
                    partners[first] = None
            else:   #first and second were both non-empty: check partners
                if first in partners:
                    if partners[first] != second:
                        self._conflicting = tuple(sorted((first, second)))
                        return True
                if second in partners:
                    if partners[second] != first:
                        self._conflicting = tuple(sorted((first, second)))
                        return True
                #add current pair to the list of constraints
                partners[first] = second
                partners[second] = first
        #can only get here if there weren't conflicts
        return False
            
    def mismatches(self, sequence, pairs=None):
        """Counts pairs that can't form in sequence.

        Sequence must have a Pairs property that acts like a dictionary
        containing a 2-element tuple for each valid pair. Can also pass in
        the pairs explicitly.
        """
        mismatches = 0
        if pairs is None:
            try:
                pairs = sequence.Alphabet.Pairs
            except AttributeError:
                pairs = sequence.Pairs
            
        for up, down in self.directed():
            curr = (sequence[up], sequence[down])
            if curr not in pairs:
                mismatches += 1
        return mismatches

    @property
    def pseudoknots(self):
        if self._pseudoknots is None:
            self._pseudoknots = self.getPseudoknots()
        
        return self._pseudoknots

    def hasPseudoknots(self):
        """Returns True if the pair list contains pseudoknots.
        
        (good_up,good_down) <=> (checked_up,checked_down)
        pseudoknot if checked_up<good_down and checked_down>good_down
        """
        return bool(self.pseudoknots)
    
    def getPseudoknots(self):
        """Returns a dict with pseudoknots. e.g. {(1, 12): [(11, 13)]}
        
        The following code was based on 
        cogent.struct.rna2d.Pairs.hasPseudoknots"""
        pseudoknots = set()
        pairs = self.directed()
        pairs.sort()
        for pos, pair in enumerate(pairs):
            if pos > 0:
                to_be_checked = pairs[:pos] + pairs[pos + 1:]
            else:
                to_be_checked = pairs[pos + 1:]
            
            for to_be_c in to_be_checked:
                if to_be_c[0] < pair[0] and pair[0] < to_be_c[1] < pair[1]:
                    pseudoknots.add(pair)
                    
        return sorted(pseudoknots)

        
    def toPartners(self, length, offset=0):
        """Returns a Partners object, if possible.
        
        length of resulting sequence must be specified.
        offset is optional, and is added to each index.
        """
        result = [None] * length
        for up_, down_ in self:
            upstream = up_ + offset
            downstream = down_ + offset
            
            if result[upstream] is not None:
                result[upstream].append(downstream)
            else:
                result[upstream] = [downstream]
                
            if result[downstream] is not None:
                result[downstream].append(upstream)
                
            else:
                result[downstream] = [upstream]
            
        result = map(lambda x: tuple(x) if x is not None else None, result)
                
        return Partners(result)
        
    def remove_conflicting_pairs(self):
        """eliminate conflicting pairs"""
        newlist = []
        visited = {}
        for a, b in self:
            if not visited.has_key(a) and not visited.has_key(b):
                newlist.append((a, b))
                visited[a] = True
                visited[b] = True
        while len(self)>0: self.pop()
        for n in newlist: self.append(n)
            
    def toVienna(self, length, offset=-1, pseudoknot_tokens=PSEUDOKNOTS_TOKENS):
        """Returns a Vienna structure string, if possible.

        length: length of resulting sequence must be specified. Instead of 
                parsing the sequence length, you can also throw in
                an object that has the required length 
                (such as the sequence that the structure corresponds to).
        offset: is optional, and is added to each index. It's set by default
                to -1 -> this is because numbering returned by RNAView starts
                from 1, not from 0 as usually in Python. And the output of
                RNAView is usually fed into BasePairs! Check calc.py.
                So for [(1,3), (4,5)] we want to have string like this: '(.)()'.
        pseudoknot_tokens: tokens used to represet pseudoknots in Vienna
                secondary structure. Default is PSEUDOKNOT_TOKENS imported
                from CONFIG.py.
        
        The following code was based on cogent.struct.rna2d.Pairs.toVienna.
        """
        if self.hasConflicts():
            raise ConflictInBasePairsError("RNA structure in Vienna format \
cannot be generated! '%s' contains conflicting base pairs!" % str(self))
        try:
            length = int(length)
        except ValueError: #raised when length can't be converted to int
            length = len(length)

        p = self.directed()
        result = ['.'] * length
        for up, down in p:
            try:
                upstream = up + offset
                downstream = down + offset
            except TypeError:
                continue
            result[upstream] = '('
            result[downstream] = ')'
        
        tokens_gen = cycle(pair for pair in pseudoknot_tokens)
        
        previous = None
        tokens = tokens_gen.next()
        for pseudoknot in self.pseudoknots:
            
            if previous is not None and not\
            (previous[0] < pseudoknot[0] < pseudoknot[1] < previous[1] \
             or pseudoknot[0] < previous[0] < previous[1] < pseudoknot[1]):
                # kolejny pseudowęzeł (niezagnieżdżone pary zasad),
                # dlatego dobieramy nowy token
                
                between = result[pseudoknot[0] : pseudoknot[1]]
                if between.count(tokens[0]) != between.count(tokens[1]):                    
                    tokens = tokens_gen.next()

            result[pseudoknot[0] + offset] = tokens[0]
            result[pseudoknot[1] + offset] = tokens[1]            

            previous = pseudoknot
            
        vienna = ViennaStructure(''.join(result))
        
        pairs_for_validation = vienna.toPairs()
        pairs_for_validation_diffs = [x[0] - x[1] for x in pairs_for_validation]
        self_diffs = [x[0] - x[1] for x in sorted(self)]
        
        if pairs_for_validation_diffs != self_diffs:
            raise PseudoknotTokensError('Did not correctly transform '+\
                                        'self -> vienna -> self '+\
                                        str(vienna) + ' ' + str(self) +\
                                        ' Probably not enough pseudoknots '+\
                                        'tokens. ' + str(pseudoknot_tokens))
        return vienna
        
    def make_non_conflicting_bps(self):
        return solve_conflicts(copy.deepcopy(self))
        
    def make_non_conflicting_viennas(self, length, offset=-1,
                                     pseudoknot_tokens=PSEUDOKNOTS_TOKENS):
        
        non_conf_conf = self.make_non_conflicting_bps()
        non_conf, conflictss = non_conf_conf[0], non_conf_conf[1:]
        
        data = [non_conf.toVienna(length, offset, pseudoknot_tokens)]
        for conflicts in conflictss:
            for conflict in conflicts:
                data.append(BasePairs([conflict]).toVienna(length, offset,
                                                           pseudoknot_tokens))
        return tuple(data)
        
class ViennaStructure(StructureString):
    """Defines what characters can be used to make RNA Vienna sec struct
    
    Author: Kristian Rother
    """
    # keys are used to initialize Alphabert dict
    keys = ''
    #dict of symbols that start base pairs
    StartSymbols = {}
    #dict of symbols that end base pairs
    EndSymbols = {}
    
    for pair in PSEUDOKNOTS_TOKENS:
        keys += pair[0] + pair[1]
        StartSymbols[pair[0]] = pair[1]
        EndSymbols[pair[1]] = pair[0]
    keys += '(.)'
    StartSymbols['('] = ')'
    EndSymbols[')'] = '('
    Alphabet = dict.fromkeys(keys)
    
def Vienna(data,Energy=None):
    """Tries to extract structure and energy from string data.

    Returns (structure, energy) where energy might be None.

    structure is just anything before the first space: doesn't validate.
    
    Author: Kristian Rother
    """
    pieces = data.strip().split(None, 1)
    if not pieces:
        return ViennaStructure('', Energy)
    else:
        if not Energy:
            try:
                energy = float_from_string(pieces[1])
            except (TypeError, ValueError, IndexError):
                energy = Energy
        else: #energy given by user overrules the one in structure
            energy = Energy
    return ViennaStructure(pieces[0], energy)
    
def solve_conflicts(bps):
    """
        Author: Tomasz Puton
    """
    potentially_conflicting_bps = copy.deepcopy(bps)
    bps.remove_conflicting_pairs()
    conflicting_bps = BasePairs(list(set(potentially_conflicting_bps) - \
                                     set(bps)))
    if conflicting_bps.hasConflicts():
        return [bps] + solve_conflicts(conflicting_bps)
    else:
        return [bps, conflicting_bps]
        
    
