#! /usr/bin/env python
#-*- coding: utf-8 -*-

"""
Parsers for RNA secondary structure formats.
"""

__author__ = "Tomasz Puton"
__credits__ = "Ewa Tkalińska, Łukasz Kozłowski, Kristian Rother"
__license__ = "GNU GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

import tempfile, os,  os.path, re, glob, copy, string
from operator import itemgetter

from cogent.parse.ct import ct_parser
from cogent.parse.stockholm import MinimalStockholmParser
from moderna.ModernaSequence import Sequence

from knotted2nested.ct_simple import MinimalCtParser
from secstruc import BasePairs, ViennaStructure, PSEUDOKNOTS_TOKENS, \
PseudoknotTokensError

"""
This modules contains parses for output formats of secondary structure prediction programs
for single sequences or groups of similar sequences (comparative prediction).
They all can be parsed into lists of basepairs,
and lists of base pairs can be converted to dot-bracket sequences

Example of structure without pseudoknots, takes one line:
    ((..(..)..).((...).(((...)))).)
Example of structure with pseudoknots, takes more than one line:
    ((..(..)..)..((...).(((...)))).)
    .......((........).)............
    ...........((..........)).......
"""


class PseudoknotError(Exception): pass

fp = os.path.split(__file__ )[:-1][0]
K2N_PATH = os.path.join(fp,'knotted2nested') + os.sep

"""

Parsers for comparative methods.

These classes are used to parse output of comparative predictions methods.
They return RNA secondary structure in vienna notation
    cons_struc - consensus structure for given sequence
    aln_strucs - for each sequence their name and structure,
             where structures are formed to be able to align with each other
    cons_level_strucs - for each sequence their name and structure,
                    where all structure represent the same abstract shapes
    differ_levels_strucs - mfe structures represent different level of abstract shapes
and sequence in form of:
    cons_seq - consensus sequence for given sequences
    aln_seqs - for each sequence their name and structure including brake-lines
"""
#######################
class ParserStockholmFormat:
    """
    Parser for output in Stockholm format from RSmatch program.
    Class return consensus structure and sequence alignment.
    Input is file or string.
    """
    def __init__(self, textfilepath):
        try:
            self.text = open(textfilepath).read()
        except:
            self.text = textfilepath
        self.cons_struc = []
        self.aln_seq = []

    def parse_stockholm_format(self):
        """
        Function cut from text part in Stockholm format and return it as a list of rows.
        Text in this format starts with '# STOCKHOLM 1.0' and ends with '//'.
        """
        Stockholm_format_text = self.text.strip().split('# STOCKHOLM 1.0')[1].split('//')[0]
        Stockholm_lines = Stockholm_format_text.strip().split('\n')
        return self.perform_parsing(Stockholm_lines)

    def perform_parsing(self, Stockholm_lines):
        """
        Function parse sequence alignment which is returned as a list and
        consensus structure from lines starts with '#=GC SS_cons' returned as a list.
        """
        index = 0
        while not Stockholm_lines[index].startswith('#=GC SS_cons'):
            self.aln_seq.append(' '.join (Stockholm_lines[index].split()) )
            index += 1
        numb_aln_seq = index
        self.cons_struc.append(Stockholm_lines[index].split()[2])
        index += 2
        while index < len(Stockholm_lines):
            for add_index in range(numb_aln_seq):
                self.aln_seq[add_index] += Stockholm_lines[index+add_index].split()[1]
            self.cons_struc[0] += Stockholm_lines[index+add_index+1].split()[2]
            index += add_index+3
        return self.aln_seq, self.cons_struc


class ParserRNAforester:
    """
    Parser for output of RNAforester multiple alignment.
    It return  alignment of sequence and structure, consensus sequence and structure
    for each cluster of sequences detected in inputfile.
    """
    def __init__(self, textfile):
        try:
            self.text = open(textfile).read()
        except:
            self.text = textfile

    def parse(self):
        """
        Function split text from file into parts represent different clusters and
        from each part retrieve alignment of sequence and structure,
        consensus sequence and structure.
        """
        result = {}
        clusters = self.text .split('RNA Structure Cluster Nr')[1:]
        for cluster in clusters:
            aln_cons = cluster.split('Consensus sequence/structure:')
            aln_seq, aln_struc = self.parse_alignment(aln_cons[0])
            cons_seq, cons_struc = self.parse_consensus(aln_cons[1])
            yield aln_seq, aln_struc, cons_seq, cons_struc

    def parse_alignment(self, alignment):
        """
        Function is used by partition_textfile to retrieve alignments.
        """
        aln_seq = []
        aln_struc = []

        groups = alignment.strip().split('\n\n')[1:]
        length = len(groups)
        aln_seq = self.parse_aln_seq(groups[:int(length/2)])
        aln_struc = self.parse_aln_struc(groups[int(length/2):])
        return aln_seq, aln_struc

    def parse_aln_seq(self, groups):
        """
        Function is used by parse_aln to retrieve alignment of sequences.
        """
        aln_seq = []
        for x in groups[0].split('\n')[:-1]:
            aln_seq.append('')
        for group in groups:
            for index, seq in enumerate(group.split('\n')[:-1]):
                aln_seq[index] += seq[24:].strip()
        return aln_seq

    def parse_aln_struc(self, groups):
        """
        Function is used by parse_aln to retrieve alignment of structures.
        """
        aln_struc = []
        for x in groups[0].split('\n')[:-1]:
            aln_struc.append('')
        for group in groups:
            for index, ss in enumerate(group.split('\n')[:-1]):
                aln_struc[index] += ss[24:].strip()
        return aln_struc

    def parse_consensus(self, consensus):
        """
        Function is used by partition_textfile to retrieve consensus sequence and structure.
        """
        groups = consensus.strip().split('\n\n')
        seq_struc = groups[0].split('\n')[10:12]
        cons_seq = seq_struc[0].strip()
        cons_struc = seq_struc[1].strip()
        for group in groups[1:]:
            seq_struc = group.split('\n')[10:12]
            cons_seq += seq_struc[0].strip()
            cons_struc += seq_struc[1].strip()
        return cons_seq, cons_struc

class ParserMASTR:
    """
    Parser for program MASTR output.
    It return  alignment of sequences and consensus structure.
    """
    def __init__(self, textfilepath):
        try:
            self.text = open(textfilepath).read()
        except:
            self.text = textfilepath
        self.aln_seq = []
        self.cons_struc = ''

    def parse_text(self):
        aln_seq = []
        lines = self.text.split('Evaluate optimized alignment:\n')[1].split('\n')
        index = 0
        while not lines[index].startswith('>structure'):
            if index%2 == 0:
                self.aln_seq.append(lines[index][1:])
            index += 1
            if index%2 == 1:
                self.aln_seq.append(lines[index])
            index += 1
        index += 1
        self.cons_struc = lines[index]
        return self.aln_seq, self.cons_struc

class ParserFoldalign:
    """
    Parser for output of Foldalign program.
    It return  alignment of two sequences and consensus structure
    for each pair comparison.
    """
    def __init__(self, textfilepath):
        try:
            self.text = open(textfilepath).read()
        except:
            self.text = textfilepath
        self.aln_seq = ['','']
        self.cons_struc = ['']

    def separate_comparison(self):
        """
        String from inputfile is splited into independent part
        describing sequences pair comparison.
        """
        list_part = self.text.strip().split('; REALIGNING')
        return list_part[1:]

    def parse_text(self, part=None):
        """
        Function parse for one part delivered by separate_comparison and for each return
        consensus sequence and sequences alignment.
        """
        if part:
            self.text = part
        text_lines = self.text.split('\n')
        for line in text_lines[:]:
            if not line.startswith('; ALIGN '):
                text_lines.remove(line)
        length = len(text_lines)
        self.aln_seq = ['', '']
        self.cons_struc = ['']
        self.aln_seq[0] += text_lines[7].split()[2]+' '
        self.aln_seq[1] += text_lines[9].split()[2]+' '
        index = 7
        split_line = text_lines[index]
        while index <= length-6:
            self.aln_seq[0] += ( ''.join(text_lines[index].split()[3:]) )
            self.cons_struc[0] += (''.join(text_lines[index+1].split()[3:]) )
            self.aln_seq[1] += ( ''.join(text_lines[index+2].split()[3:]) )
            index+= 4
        return self.aln_seq, self.cons_struc

class ParserRNAcast:
    """Parser for output of RNAcast program.
    It return  sequences and their names and structures,
    where all structures are the same from the chosen level of abstract.
    """
    def __init__(self, textfilepath):
        try:
            self.text = open(textfilepath).read()
        except:
            self.text = textfilepath
        self.cons_level_strucs = []

    def parse_text(self):
        lines_list = self.text.strip().split('\n')[1:]
        index = 0
        while index < len(lines_list):
            name = lines_list[index][1:]
            seq = lines_list[index+1].strip()
            struc = lines_list[index+2].split()[1]
            self.cons_level_strucs.append(name+' '+seq+' '+struc)
            index += 3
        return self.cons_level_strucs

##                                                                           ##
## Parsers for single secondary structure predictions.                       ##
##                                                                           ##

def parse_bp(data):
    """Returns a list of base pairs parsed from a BP format line iterator."""
    
    header = []
    out_data = []
    seq = ""
    
    if type(data) in (str, unicode):
        data = data.split('\n')
    elif type(data) is file:
        pass
    elif type(data) in (tuple, list):
        pass
    else:
        assert TypeError('input is neither a file nor a string/unicode!')
    
    for line in data:
        
        line = line.strip()
        if len(line) == 0:
            continue
        if not line[0] in string.digits:
            header.append(line)
        else:
            items = line.strip().split()
            assert len(items) == 3, items
            up_, down_ = (int(items[0]), int(items[2]))
            seq += items[1]
            if down_ == 0:
                continue
            
            sorted_pair = tuple(sorted((down_, up_)))
            if sorted_pair not in out_data:
                out_data.append(sorted_pair)
            
    out_data.sort()
    return BasePairs(out_data)

def parse_rnashapes(fh):
    """Returns a list of base pairs parsed from a RNAshapes output format
        line iterator.
        
        fh: file handle or anything behaving like it
        """
    data = [d for d in fh]
    start_index = 1
    if data[0].startswith('>'):
        start_index += 1

    for line in data[start_index:]:
        tokens = line.strip().split()
        if len(tokens)==3:
            energy, secstruc, shape = tokens
            yield float(energy), secstruc, shape

def parse_hotknots_v20(open_file, seq_len):
    """A parser for HotKnots in version 2.0 and higher
    """
    for line in open_file:
        if line.startswith('Total number of RNA structures:'):
            seq_line = open_file.next()
            while True:
                ss_line = open_file.next().split()
                vienna = ss_line[1]
                energy = float(ss_line[2])
                yield (vienna, energy)

def parse_hotknots(open_file, seq_len):
    """A very non elegant solution to parsing HotKnots output.
    Alternative structs are completely ignored.
    """
    for line in open_file:
        if line.startswith(' total number of Rna structures:'):
            while True:                
                energy_line = open_file.next()
                try:
                    energy = float(energy_line.split()[-2])
                except IndexError:
                    raise StopIteration
                
                sec_struc_line = open_file.next()
                vienna = '.' * seq_len
                
                pseudoknots_tokens = iter([('[',']'), ('{','}'), ('<','>')])
                
                for block in sec_struc_line.split(',')[:-1]:
                    start_range, end_range = \
                        [x.strip() for x in block.split(';')]
                    split_start = start_range.split('-')
                    split_end   = end_range.split('-')
                    
                    if not len(vienna) == vienna.count('.') and vienna[int(split_start[1]) : int(split_end[0]) + 1].count(')') == vienna[:int(split_start[0]) + 1].count('('):
                        tokens = pseudoknots_tokens.next()
                    else:
                        tokens = ['(',')']
                    
                    vienna_start_part = vienna[:int(split_start[0])-1] + tokens[0] * len(range(int(split_start[0]), int(split_start[1])+1)) + vienna[int(split_start[1]):]
                    vienna = vienna_start_part + vienna[len(vienna_start_part):]
                    vienna_end_part = vienna[:int(split_end[0])-1] + tokens[1] * len(range(int(split_start[0]), int(split_start[1])+1)) + vienna[int(split_end[1]):]
                    vienna = vienna_end_part
                    
                yield (vienna, energy)

def parse_pknots(fh):
    """Returns a list of base pairs and their corresponding energy parsed
        from a pknots format line iterator.
        
        fh: a file open for reading
    """
    energy = None
    last_number = '-1'
    pairs_str = []
    for line in fh:
        stripped = line.strip()
        if stripped.startswith('0') or \
        stripped.startswith(str(int(last_number) + 1)):
            
            first_line = stripped.split()
            second_line = fh.next().strip().split()
            for op, cl in zip(first_line, second_line):
                if cl == '.':
                    continue
                if (cl, op) not in pairs_str:
                    pairs_str.append((op, cl))
            last_number = first_line[-1]
        
        if stripped.startswith('energy'):
            energy = float(line.split()[-1].strip())
            break
        
    return (BasePairs([(int(x)+1, int(y)+1) for x,y in pairs_str]), energy)

def parse_vienna(data, row_by_row=False):
    """
    Returns a secondary structure in the row indicated taking care of getting
    energy also.
    """
    data = [d for d in data]
    energy_regex = re.compile(r'.\d+\.\d*')
    if row_by_row:
        seq_name = None
        seq = None
        ss = data[0].split()[0].strip()
        res = energy_regex.search(data[0])
    else:
        seq_name = data[0][1:]
        seq = data[1].strip()
        ss = data[2].split()[0].strip()
        res = energy_regex.search(data[2])
        
    if res:
        tmp = res.group().replace('(','').strip()
        energy = float(tmp)
    else:
        energy = None
        
    return (seq_name, seq, ss, energy)

def parse_multiple_vienna(data, row_by_row=True):
    """
    Yields tuples of dot-bracket strings and energies from a multiline vienna
    format line iterator.
    """
    data = [d for d in data[2:]]
    for line in data:
        yield parse_vienna([line], row_by_row=row_by_row)

def parse_ct(data):
    """Returns a list of base pairs parsed from a CT format line iterator."""
    pairs = [r[2] for r in MinimalCtParser(data)][0]
    result = [(x+1, y+1) for x, y in pairs]
    result.sort()
    return BasePairs(result)

def parse_afold(data):
    """Returns a tuple consisting of a list of base pairs parsed from an Afold
    format line iterator and an energy for that folding."""
    bpairs = []
    energy = None
    for line in data:
        if line.find('Multidomain') != -1:
            energy = float(line.split('=')[1].split()[0].strip())
        if line.strip()=='': continue
        first = line.split()[0]
        tokens = line.replace('.', ' ').split()
        bpairs.append((int(tokens[1]), int(tokens[2])))
    return (BasePairs(bpairs), energy)

def parse_sfold(data):
    """Returns a list of base pairs parsed from an Sfold format line iterator."""
    bpairs = []
    for line in data:
        if re.search("\d+\s+\d+\n*\Z", line):
            a, b = line.strip().split()
            bpairs.append((int(a), int(b)))
    return BasePairs(bpairs)

def is_bp_a_pseudoknot(bp0, bp1, bplist):
    """Returns True if a given base pair is a pseudoknot, False otherwise.

    bp0: int, number of the first base
    bp1: int, number of the second base
    bplist: a list of base pair indices
    """
    for pair in BasePairs(bplist).directed():
        if bp0 < pair[0] < bp1 < pair[1] or pair[0] < bp0 < pair[1] < bp1:
            return True
    return False

def make_dotbracket_from_bplist(length, bplist, \
    pseudoknot_tokens=PSEUDOKNOTS_TOKENS):
    """
    Takes a structure as a list of basepairs [(2,7),(4,5),(10,15),(17,20)]
     and converts it to vienna format: "((...).(((.((....))).)))"

    length: int, length of an RNA sequence

    base_pairs: a list of lists - every sub-list contains two integers which are
    base numbers (numbering starts from 1!!! not from 0)

    pseudoknot_tokens: a list. Pseudoknot tokens are used to represent
    pseudoknots in RNA secondary structure in a Vienna format. The user is
    responsible for providing sufficient number of pseudoknot tokens. Usually
    [('[',']'), ('{','}'), ('<','>')] will be enough, but in situations when
    an RNA contains more than 3 pseudoknots, using the default
    pseudoknot_tokens will result in PseudoknotTokensError exception.
    """
    base_pairs = BasePairs(bplist)
    return base_pairs.toVienna\
        (length=length, pseudoknot_tokens=pseudoknot_tokens, offset=-1)

def make_bplist_from_dotbracket(secstruc):
    """Reads a dot-bracket secondary structure and returns BasePairs object.

    secstruc: str, a dot-bracket secondary structure.
    """
    vienna = ViennaStructure(secstruc)
    return vienna.toPairs()

def write_bp_file(bplist, filename):
    """Dirty hack bp file writer."""
    # prepare base pair list for writing bp file
    bplist = bplist[:]
    new = []
    for bp in bplist: new.append((bp[1], bp[0]))
    bplist += new
    maxi = max([bp[1] for bp in bplist])
    i = 1
    while i<maxi:
        found = False
        for bp in bplist:
            if bp[0]==i: found = True
        if not found:
            bplist.append((i, 0))
        i += 1
    bplist.sort()
    # write bp file
    bpdata = ["Filename: example_small.bpseq\n",
    "Organism: Some organism\n",
    "Accession Number: XYZ123\n",
    "Citation and related information available at http://www.rna.ccbb.utexas.edu\n"]
    bpdata += ["%i A %i\n"%(bp[0], bp[1]) for bp in bplist]
    open(filename, 'w').writelines(bpdata)

def has_overlapping_basepairs(bplist):
    indices = {}
    for bp in bplist:
        if indices.has_key(bp[0]) or indices.has_key(bp[1]): return True
        indices[bp[0]]=True
        indices[bp[1]]=True
    return False

def remove_pseudoknot(bplist):
    """Returns a list of base pairs, with pseudoknots removed but the
    highest possible number of base pairs preserved."""
    if has_overlapping_basepairs(bplist):
        raise PseudoknotError('Base pairs have duplicate indices. Cannot remove pseudoknot!')
    result = []
    file_id_in,  in_bp_name = tempfile.mkstemp(suffix='.bp', prefix='pkremoval_in_')
    write_bp_file(bplist, in_bp_name)
    file_id_out,  out_bp_name = tempfile.mkstemp(suffix='.bp', prefix='pkremoval_out_')
    os.system('python %sknotted2nested.py -f bpseq -F bpseq %s > %s'%(K2N_PATH, in_bp_name,  out_bp_name))
    result = parse_bp(open(out_bp_name))
    os.close(file_id_in)
    os.unlink(in_bp_name)
    os.close(file_id_out)
    os.unlink(out_bp_name)


    return result

class RNAStrandRecordException(Exception): pass

class RNAStrandRecord(object):
    """Class representing RNAstrand CT record.
    """
    def __init__(self, ct_path):
        """
        ct_path: str, a path to RNAStrand's ct file
        """
        self.ct_path = ct_path
        self._defline, self._sequence, self._structure, self._pairs = \
            "", "", "", ""
    
    def __repr__(self):
        return '<%s RNA STRAND record>' % self.ct_path
    
    @staticmethod
    def _is_valid_sequence(raw_sequence):
        """Validates sequence extracted from RNAstrand data file in
        dot-parentheses format.
        """
        if re.search(r'[_XN~?]', raw_sequence) or len(raw_sequence) <= 10:
            return False
        return True
    
    def _parse(self):

        tmp_ct = []
        defline = ">"
        for line in open(self.ct_path):
            if line.startswith('#'):
                defline += line.strip() + '|'
            else:
                tmp_ct.append(line)
        defline = defline.replace('\n','').replace('#','')
        
        # just one entry per file
        result = ct_parser(tmp_ct)[0]
        sequence = result[0]
        
        # IMPORTANT!
        # ct_parser starts numbering of base pairs from 0, not from 1, as
        # elsewhere in nucleic.secstruc. To make things consistent, here
        # we'll add 1 to each index in base pairs list!
        result_1 = [(pair[0]+1, pair[1]+1) for pair in result[1]]
        pairs = BasePairs(result_1)
        
        if self._is_valid_sequence(sequence):
            sequence = Sequence(sequence.upper()).seq_without_modifications
            # Second check, because ModeRNA might also mix a bit at this stage
            # including for exampel X in sequences without modifications!
            if self._is_valid_sequence(sequence):
                try:
                    vienna = pairs.toVienna(len(sequence))
                except PseudoknotTokensError:
                    vienna = None
                finally:
                    return defline, sequence, vienna, pairs
            
        return None, None, None, None
    
    def _parse_if_necessary(self):
        if not \
            (self._defline or self._sequence or self._structure or self._pairs):
            self._defline, self._sequence, self._structure, self._pairs \
                = self._parse()
        
    @property
    def valid(self):
        self._parse_if_necessary()
        if self.sequence is None:
            return False
        else:
            return True
        
    @property    
    def defline(self):
        self._parse_if_necessary()
        return self._defline
    
    @property
    def pairs(self):
        self._parse_if_necessary()
        return self._pairs
        
    @property
    def sequence(self):
        self._parse_if_necessary()
        return self._sequence

    @property
    def structure(self):
        self._parse_if_necessary()
        return self._structure
        
    @property
    def vienna(self):
        """A Vienna representation of RNAstrand record.
        Defline is created using values found after #.
        """
        self._parse_if_necessary()
        if self._sequence and self._structure and self._defline:
            return self._defline + '\n' + self._sequence + '\n' + \
            self._structure +'\n'
    
    @property    
    def fasta(self):
        """A FASTA representation of RNA sequence from RNAstrand record.
        Defline is created using values found after #.
        """
        self._parse_if_necessary()
        if self._sequence and self._structure and self._defline:
            return self._defline + '\n' + self._sequence + '\n'
    
class RNAStrandParser(object):
    """Parser for a big file with RNAstrand records"""
    
    def __init__(self, path, sorted_by_size=False):
        """path: str, a path to a directory with RNAstrand records in 
        files with modified CT format.
        """
        self.path = path
        self.ct_files = glob.iglob('%s/*.ct' % self.path)
        if sorted_by_size:
            self.ct_files = iter(sorted(self.ct_files, cmp=self.__size_sorter__))
            
    @staticmethod
    def __size_sorter__(x, y):
        return cmp(os.path.getsize(x), os.path.getsize(y))
    
    def __repr__(self):
        return 'RNA STRAND parser in %s' % self.path
    
    def __iter__(self):
        return self
    
    def next(self):
        return RNAStrandRecord(self.ct_files.next())
        
    def make_length_histogram(self):
        """Makes length dependent histogram using matplotlib
        """
        seqs_lengths = []
        for rec in self:
            if rec.sequence:
                seqs_lengths.append(len(rec.sequence))
        bins = range(1, 4000, 10)
        try:
            from matplotlib import pylab
        except ImportError:
            print "This feature requires having matplotlib installed!"
        else:
            pylab.hist(seqs_lengths, bins)
            pylab.show()
        
     
class Token(object):
    def __init__(self, token_open, token_close):
        self._token_open = str(token_open)
        self._token_close = str(token_close)
        
    def __str__(self):
        return self._token_open
    
    def __repr__(self):
        return repr(self._token_open)
        
    def __eq__(self, other):
        return self.opening == other
    
    @property
    def opening(self):
        return self._token_open
    
    @property
    def closing(self):
        return self._token_close
     
# letters from a to z
UNACCEPTED_TOKENS = [Token(op, cl) for (op, cl) in \
                     [(chr(i), chr(i).upper())for i in range(97, 123)]]

# letters from A to Z
UNACCEPTED_TOKENS += [Token(op, cl) for (op, cl) in \
                     [(chr(i), chr(i).upper())for i in range(65, 91)]]

# < and > ;)
UNACCEPTED_TOKENS += [Token('<', '>')]

GOOD_TOKENS = [Token(op, cl) for (op, cl) in [('{','}'), ('[',']')]] * 100
        
        
def fix_pseudoknot(bad_tokens, good_tokens, vienna, block_len):
    """Fixes and returns RNA secondary structure string containing pseudoknots.
        
        Takes a 'vienna' RNA secondary structure string and replaces
        'bad_tokens' with 'good_tokens'.
        Each 'block_len' occurances of bad token will be replaced.
        
        bad_tokens: a tuple of two one character strings corresponding to an
            opening and closing token e.g. ('a', 'A')
        good_tokens: a tuple of two one character strings corresponding to an
            opening and closing token e.g. ('[', ']')
        vienna: str, an RNA secondary structure in Vienna format with
            pseudoknot tokens to be fixed
        block_len: int, a 4th parameter passed to re.sub (see re.sub docs);
            it sets the number of token occurences to be replaced
            (setting this parameter correctly is crucial for this stuff to work)
            Check the doctests below!
        
        >>> fix_pseudoknot(Token('a','A'), Token('[', ']'), 'aa..AA', 2)
        '[[..]]'

        >>> fix_pseudoknot(Token('A','A'), Token('[', ']'), 'AA..AA', 2)
        '[[..]]'
        
        >>> fix_pseudoknot(Token('A','A'), Token('[', ']'), 'AB..AB', 1)
        '[B..]B'
        
        >>> fix_pseudoknot(Token('B','B'), Token('(', ')'), '[B..]B', 1)
        '[(..])'
        
        Returns a vienna string with replaced bad tokens into good tokens.
    """
    return re.sub(bad_tokens.closing, good_tokens.closing, \
                  re.sub(bad_tokens.opening, good_tokens.opening, \
                         vienna, block_len),\
                  block_len)
    
def identify_bad_tokens(vienna):
    """
    Identifies bad pseudoknot tokens and returns a list of tuples.
    Each tuple contains two characters (one corresponding to an
    opening token and the other to the closing token; they can be the same).
    
    vienna: str, an RNA secondary structure in Vienna format with
        pseudoknot tokens to be fixed
    
    See the doctest below.
    
    >>> identify_bad_tokens('.....')
    []
    
    >>> identify_bad_tokens('(()).')
    []

    >>> identify_bad_tokens('A...A')
    ['A']

    >>> identify_bad_tokens('a().A')
    ['a']

    >>> identify_bad_tokens('aaaaa().AAAAA')
    ['a']

    >>> identify_bad_tokens('aaaaa(bbb).AAAAABBB')
    ['a', 'b']

    >>> identify_bad_tokens('A.(.A.).')
    ['A']
    
    >>> identify_bad_tokens('AAAA.(.AAAA.).')
    ['A']
    
    >>> identify_bad_tokens('AAAABBBB.(.AAAABBBB.).')
    ['A', 'B']
    
    >>> identify_bad_tokens\
        ('...AAA.(((((.........AAA.......)))))..((((.....))))..')
    ['A']
    
    """
    skipped = []
    bad = []
    for i in vienna:
        if i in UNACCEPTED_TOKENS and i not in skipped:
            u_token = UNACCEPTED_TOKENS[UNACCEPTED_TOKENS.index(i)]
            bad.append(u_token)
            skipped.append(u_token.opening)
            skipped.append(u_token.closing)
    return bad

def has_weired_pseudoknots(vienna):
    """Checks whether a vienna string has bad pseudoknots' characters.
    
        vienna: str, RNA secondary structure in Vienna format e.g. ..(..)..
    
        >>> has_weired_pseudoknots('.()..')
        False
        >>> has_weired_pseudoknots('.A(A).')
        True
        >>> has_weired_pseudoknots('.a(A).')
        True
        >>> has_weired_pseudoknots('.aaa(AAA).')
        True
        >>> has_weired_pseudoknots('..<.(..>).')
        True
    """
    for i in vienna:
        if i in UNACCEPTED_TOKENS:
            return True
    return False
    
class PseudoknotFormattingError(Exception): pass
        
def unify_pseudoknot_chars(vienna, good_tokens=GOOD_TOKENS, bad_tokens=[]):
    """Unifies pseudoknot characters in Vienna string
    
        vienna: str, RNA secondary structure in Vienna format e.g. ..(..)..
        
        >>> unify_pseudoknot_chars('....', GOOD_TOKENS, [])
        '....'
        
        >>> unify_pseudoknot_chars('.().', GOOD_TOKENS, [])
        '.().'
        
        >>> unify_pseudoknot_chars('[(])', GOOD_TOKENS, \
        [Token('a', 'A'), Token('b', 'B')])
        '[(])'

        >>> unify_pseudoknot_chars('a(A)', GOOD_TOKENS, [Token('a', 'A')])
        '[(])'
        
        >>> unify_pseudoknot_chars('.aa..(.AA.).b..(.B..)', \
        GOOD_TOKENS, [Token('a', 'A'), Token('b', 'B')])
        '.{{..(.}}.).[..(.]..)'
        
        >>> unify_pseudoknot_chars\
        ('...AAA.(((((.........AAA.......)))))..((((.....))))..', \
         GOOD_TOKENS, [Token('A', 'A')])
        '...[[[.(((((.........]]].......)))))..((((.....))))..'
    """
    if len(good_tokens) < len(bad_tokens):
        raise PseudoknotFormattingError(\
            "Not enough tokens to reformat pseudoknots in '%s'. \
Good tokens: %s, bad tokens: %s." % (vienna, good_tokens, bad_tokens))
    
    good_tokens_cp = copy.deepcopy(good_tokens)
    bad_tokens_cp = copy.deepcopy(bad_tokens)
    
    if not bad_tokens:
        bad_tokens_cp = identify_bad_tokens(vienna)
    
    if bad_tokens_cp == []:
        return vienna
    
    good = good_tokens_cp.pop()
    bad = bad_tokens_cp.pop()
    if bad.opening == bad.closing:
        count = vienna.count(bad.opening) / 2
    else:
        if vienna.count(bad.opening) != vienna.count(bad.closing):
            raise PseudoknotFormattingError\
                ('Your Vienna record is a complete crap! \'%s\'' % vienna)
        count = vienna.count(bad.opening)

    vienna = fix_pseudoknot(bad, good, vienna, count)
    
    if has_weired_pseudoknots(vienna):
        return unify_pseudoknot_chars(vienna, good_tokens_cp, bad_tokens_cp)
    else:
        return vienna
    
    
    
    
def compose_ss_from_cmfinder_motives(motives, seq_name, sequence):
    """A parser for CMfinder motifs. Tries to create a secondary structure 
        
        motives: list of handles to file-like objects with a motif
        seq_name: str, name of sequence for which motives are to be extracted
    """
    assert type(seq_name) is str
    gs_regex = re.compile(r'^#=GS\s+%s\s+DE' % re.escape(seq_name))
    gr_regex = re.compile(r'^#=GR\s+%s\s+SS' % re.escape(seq_name))
    data = []

    for motif_fh in motives:
        assert hasattr(motif_fh, 'read')
        
        read = motif_fh.readlines()
        parser = MinimalStockholmParser(read)
        record = parser.next()[0]
        for gs in record['GS']:
            gs_search = gs_regex.search(gs)
            if gs_search:
                split = gs.split()
                if len(split) == 5:
                    header, seq_name_, de, start_stop, score = split
                    start, dots, stop = start_stop.split('.')
                    del dots, de
                elif len(split) == 6:
                    header, seq_name_, de, start, stop, score = split
                    start = start.split('.')[0]
                else:
                    raise ValueError('non-std output of CMfinder! %s' % \
                                     str(split))
                
                del header
                start = int(start)
                stop = int(stop)
                score = float(score)
                
                # aligned seq read from CMfinder output in Stockholm format
                aligned_seq = "".join(x.split()[1].strip() for x in read if \
                               x.startswith(seq_name))
                
                ss = ''
                for rec in record['GR']:
                    gr_search = gr_regex.search(rec)
                    if gr_search:
                        ss += rec.split()[-1]
                        
                degapped_ss  = ''
                for c_, e_ in zip(aligned_seq, ss):
                    if c_ in ['.', '-', '_']:
                        continue
                    degapped_ss  += e_
                        
                ss = degapped_ss.replace('<', '(').replace('>', ')')\
                    .replace('-','.')
                
                assert seq_name == seq_name_
                data.append((start, stop, score, ss))
            
    # 1. first sort by score and get motif with the highest score
    data.sort(key=itemgetter(2), reverse=True)
    highest_score = data[0][2]
    with_highest_scores = [(x[1] - x[0], x[0], x[1], x[2], x[3]) for x in \
        data if x[2] == highest_score]
    
    # 2. get the longest motif with that score
    with_highest_scores.sort(key=itemgetter(0), reverse=True)
    
    final_motif = with_highest_scores[0][1:]
    final_motif_coords = range(final_motif[0], final_motif[1]+1)
    
    # 3. let's see whether it is possible to add more motives to the best one,
    # provided that they don't overlap
    final_motives = [final_motif]
    for d in data:
        if d == final_motif:
            continue
        
        # only motives located before or after are allowed
        if d[1] < min(final_motif_coords) or d[0] > max(final_motif_coords):
            # ok, merging!
            final_motif_coords += range(d[0], d[1]+1)
            final_motives.append(d)
    
    final_ss = ['.'] * len(sequence)
    for motif in final_motives:
        start = motif[0]
        stop = motif[1]
        for index, ss_element in enumerate(motif[3]):
            final_ss[index+start-1] = ss_element
        
    return "".join(final_ss)
    