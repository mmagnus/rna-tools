#!/usr/bin/env python
"""Code for handling RNA secondary structure.

RNA secondary structures can be represented in many different ways. The
representation makes a large difference to the efficiency of different
algorithms, so several different structural representations (and the means
to interconvert them) are provided.

Provides following classes:  
    Stem: representation of a stem in a secondary structure
    Partners: list holding partner of each position.
    PseudoknowRemover: detects and removes pseudoknots.
    Pairs: list of base pairs in a structure.
    StructureString: string holding a secondary structure representation.
    ViennaStructure: representation of Vienna-format RNA structure.
    WussStructure: Wash U secondary structure format, handles pseudoknots.
    StructureNode: for tree representation of nested RNA structure.
"""
from string import maketrans

class keep_chars(object):
    """Returns a filter object o(s): call to return a filtered string.

    Specifically, strips out everything in s that is not in keep.
    This filter is case sensitive by default.
    """
    allchars = maketrans('','')
    def __init__(self, keep, case_sens=True):
        """Returns a new keep_chars object, based on string keep"""
        if not case_sens:
            low = keep.lower()
            up = keep.upper()
            keep = low + up
        self.delchars = ''.join([c for c in self.allchars if c not in keep])
    
    def __call__(self, s, a=None, d=None):
        """f(s) -> s, translates using self.allchars and self.delchars"""
        if a is None: a = self.allchars
        if d is None: d = self.delchars
        return s.translate(a,d)

maybe_number = keep_chars('0123456789.+-eE')

def float_from_string(data):
    """Extracts a floating point number from string in data, if possible."""
    return float(maybe_number(data))


class PairError(ValueError):
    """Base class for errors in pairing."""
    pass

class Partners(list):
    """Holds list p such that p[i] is the index of the partner of i, or None.
    
    Primarily useful for testing whether a specified base is paired and, if so,
    extracting its partner.

    Each base may have precisely 0 or 1 partners. If A pairs with B, B must
    pair with A.
    All inconsistencies will be removed by setting previous partners to None.
    Checking for conflicts and raising errors should be done in method that 
    constructs the Partners.

    If constructing by hand, should initialize with list of [None] * seq_length.
    Typically, Partners will be constructed by code from some other data. Use
    the EmptyPartners(n) factory function to get an empty Partners list of
    length n.
    """
    def __setitem__(self, index, item):
        """Sets self[index] to item, enforcing integrity constraints."""
        if index == item:
            raise ValueError, "Cannot set base %s to pair with itself." % item
        #if item already paired, raise Error or make partner unpaired
        if item and self[item]:
            self[self[item]] = None
        #if already paired, need to make partner unpaired
        curr_partner = self[index]
        if curr_partner is not None:
            list.__setitem__(self, curr_partner, None)
        #set self[index] to item    
        list.__setitem__(self, index, item)
        #if item is not None, set self[item] to index
        if item is not None:
            list.__setitem__(self, item, index)
                
    def toPairs(self):
        """Converts the partners to sorted list of pairs."""
        result = Pairs()
        for first, second in enumerate(self):
            if first < second:
                result.append((first, second))
        return result

    def _not_implemented(self, *args, **kwargs):
        """Raises NotImplementedError for 'naughty' methods.
        
        Not allowed any methods that insert/remove items or that change the
        order of the items, including things like sort or reverse.
        """
        raise NotImplementedError
    
    __delitem__ = __delslice__ = __iadd__ = __imul__ = __setslice__ = append \
    = extend = insert = pop = remove = reverse = sort = _not_implemented

def EmptyPartners(length):
    """Returns empty list of Partners with specified length."""
    return Partners([None] * length)


class Pairs(list):
    """Holds list of base pairs, each of which is a 2-element sequence.
    
    This is a very lightweight object for storing base pairs, and does not
    perform any validation. Useful as an intermediate in many different
    calculations.
    """
    def toPartners(self, length, offset=0, strict=True):
        """Returns a Partners object, if possible.
        
        length of resulting sequence must be specified.
        offset is optional, and is added to each index.
        strict specifies whether collisions cause fatal errors. if not strict
        conflicts will be removed by the Partners object.
        """
        result = EmptyPartners(length)
        for up, down in self:
            upstream = up + offset
            downstream = down + offset
            
            if result[upstream] or result[downstream]:
                if strict:
                    raise ValueError, "Pairs contain conflicting partners: %s"\
                        % self
            result[upstream] = downstream
        return result
            
    def toVienna(self, length, offset=0, strict=True):
        """Returns a Vienna structure string, if possible.

        length of resulting sequence must be specified.
            Instead of parsing the sequence length, you can also throw in
            an object that has the required length (such as the sequence that
            the structure corresponds to).
        offset is optional, and is added to each index.
        strict specifies whether collisions cause fatal errors.
        """
        if self.hasPseudoknots():
            raise Exception, "Pairs contains pseudoknots %s"%(self)
        
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

            if strict:
                if (result[upstream] != '.') or (result[downstream] != '.'):
                    raise ValueError, "Pairs contain conflicting partners: %s"\
                        % self
            result[upstream] = '('
            result[downstream] = ')'
        return ViennaStructure(''.join(result))

    def tuples(self):
        """Converts all pairs in self to tuples, in place.

        Useful for constructing dicts and for sorting (otherwise, pairs of
        different types, e.g. lists and tuples, will sort according to type
        rather than to position).
        """
        self[:] = map(tuple, self)

    def unique(self):
        """Returns copy of self omitting duplicate pairs, preserving order.

        Keeps the first occurrence of each pair.
        """
        seen = {}
        result = []
        for p in map(tuple, self):
            if p not in seen:
                seen[p] = True
                result.append(p)
        return Pairs(result)

    def directed(self):
        """Returns copy of self where all pairs are (upstream, downstream).
        
        Omits any unpaired bases and any duplicates. Result is in arbitrary 
        order.
        """
        seen = {}
        for up, down in self:
            if (up is None) or (down is None):
                continue    #omit unpaired bases
            if up > down:
                up, down = down, up
            seen[(up, down)] = True
        result = seen.keys()
        return Pairs(result)

    def symmetric(self):
        """Retruns copy of self where  each up, down pair has a down, up pair.
        
        Result is in arbitrary order. Double pairs and pairs containing None 
        are left out.
        """
        result = self.directed()
        result.extend([(down, up) for up, down in result])
        return Pairs(result)
     
    def paired(self):
        """Returns copy of self omitting items where a 'partner' is None."""
        return Pairs(filter(not_none, self))
        
    def hasPseudoknots(self):
        """Returns True if the pair list contains pseudoknots.
        
        (good_up,good_down) <=> (checked_up,checked_down)
        pseudoknot if checked_up<good_down and checked_down>good_down
        """
        pairs = self.directed()
        seen = [] # list of pairs against which you compare each time
        pairs.sort()
        for pair in pairs:
            if not seen:
                seen.append(pair)
            else:
                lastseen_up, lastseen_down = seen[-1]
                while pair[0] > lastseen_down:
                    seen.pop()
                    if not seen:
                        break
                    else:
                        lastseen_up,lastseen_down = seen[-1]
                if not seen:
                    seen.append(pair)
                    continue
                if pair[1]>lastseen_down:
                    #pseudoknot found
                    return True
                else:
                    #good pair
                    seen.append(pair)
        return False


    def hasConflicts(self):
        """Returns True if the pair list contains conflicts.

        Conflict occurs when a base has two different partners, or is asserted
        to be both paired and unpaired.
        """
        partners = {}
        for first, second in self:
            if first is None:
                if second is None:
                    continue    #no pairing info
                else:
                    first, second = second, first   #swap order so None is 2nd
            if second is None: #check first isn't paired
                if partners.get(first, None) is not None:
                    return True
                else:
                    partners[first] = None
            else:   #first and second were both non-empty: check partners
                if first in partners:
                    if partners[first] != second:
                        return True
                if second in partners:
                    if partners[second] != first:
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
    
        

class StructureString(str):
    """Base class for ViennaStructure and WussStructure. Immutable.
    
    StructureString holds a structure and a energy. By default energy is 
    set to None. 
    If you compare two StructureStrings the structure is the only important
    thing, since the energy is relative to the associated sequence.
    """
    Alphabet=None
    StartSymbols = ''      #dict of symbols that start base pairs
    EndSymbols = ''        #dict of symbols that end base pairs
    
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
        num_bases = len(self) #number of bases
        result = [None] * len(self) #array of None, one for each base
        stack = []
        start = self.StartSymbols
        end = self.EndSymbols
        for i, symbol in enumerate(self):
           if symbol in start:       #open a pair
               stack.append(i)
           elif symbol in end:     #close a pair
               curr = stack.pop()  #return and delete last element
               result[i] = curr #make i pair with the last element...
               result[curr] = i #...and the last element pair with i
               
        #test whether there are any open pairs left unaccounted for        
        if stack:
           raise IndexError, \
           "Too many open pairs in structure:\n%s" % self
        return Partners(result)

    def toPairs(self):
        """Makes list of (upstream,downstream) partners.
        
        Note that the numbering starts at 0 for the first position.
        Key will always be smaller than value.

        Result is in arbitrary order.
        """
        result = {}
        stack = []
        start = self.StartSymbols
        end = self.EndSymbols
        for i, symbol in enumerate(self):
           if symbol in start:       #open a pair
               stack.append(i)
           elif symbol in end:     #close a pair
               result[stack.pop()] = i
        #test whether there are any open pairs left unaccounted for        
        if stack:
           raise IndexError, \
           "Too many open pairs in structure:\n%s" % self
        return Pairs([(key,result[key]) for key in result])


class ViennaStructure(StructureString):
    """Contains a Vienna dot-bracket structure, possibly with energy."""
    Alphabet = dict.fromkeys('(.)')
    StartSymbols = {'(':None}      #dict of symbols that start base pairs
    EndSymbols =   {')':None}      #dict of symbols that end base pairs
 
def Vienna(data,Energy=None):
    """Tries to extract structure and energy from string data.

    Returns (structure, energy) where energy might be None.

    structure is just anything before the first space: doesn't validate.
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


