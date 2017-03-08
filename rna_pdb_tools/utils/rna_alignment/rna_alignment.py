#!/usr/bin/env python

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class RNAalignment():
    """RNA alignment

    :var self.io: ``AlignIO.read(fn, "stockholm")``

    Read more: http://biopython.org/DIST/docs/api/Bio.AlignIO.StockholmIO-module.html
    """
    def __init__(self,fn):
        self.fn = fn
        self.lines = open(fn).read().split('\n')
        self.io = AlignIO.read(fn, "stockholm")
        self.ss_cons = self.get_ss_cons()
        self.copy_ss_cons_to_all()

    def copy_ss_cons_to_all(self):
        for s in self.io:
            s.letter_annotations['secondary_structure'] = self.ss_cons
            s.ss = self.ss_cons
            s.ss_clean = self.get_clean_ss(s.ss)
            s.seq_nogaps = str(s.seq).replace('-', '')
            s.ss_nogaps = self.get_ss_remove_gaps(s.seq, s.ss_clean)

    def copy_ss_cons_to_all_editing_sequence(self,seq_id, before, after):
	"""Change a sequence's sec structure.
	
	:param seq_id: string, sequence id to change, eg: 'AE009948.1/1094322-1094400'
	:param before: string, character to change from, eg: `,`
	:param after: string, character to change to, eg: `.`
	.. warning:: before and after has to be one character long
        """
        for s in self.io:
	    if s.id == seq_id
	      s.letter_annotations['secondary_structure'] = self.ss_cons.replace(before, after)
	    else:
	      s.letter_annotations['secondary_structure'] = self.ss_cons

    def get_ss_remove_gaps(self, seq, ss):
        """
        :param seq: string, sequence
        :param ss: string, ss

        UAU-AACAUAUAAUUUUGACAAUAUGG-GUCAUAA-GUUUCUACCGGAAUACC--GUAAAUAUUCU---GACUAUG-UAUA-
        (((.(.((((,,,(((((((_______.))))))).,,,,,,,,(((((((__.._____))))))...),,)))).)))).
        """
        nss = ''
        for i,j in zip(seq,ss):
            if i != '-':
                nss += j
        return nss

    def get_ss_cons(self):
        for l in open(self.fn):
            if l.startswith('#=GC SS_cons'):
                return l.replace('#=GC SS_cons','').strip()
        
    def get_shift_seq_in_align(self):
        """RF_cons vs '#=GC RF' ???"""
        for l in open(self.fn):
            if l.startswith('#=GC RF'):
                # #=GC RF                        .g.gc.a
                l = l.replace('#=GC RF','')
                c = 12 # len of '#=GC RF'
                #                         .g.gc.a
                for i in l:
                    #print i
                    if i == ' ':
                        c += 1 
                return c

    def map_seq_on_seq(self, seq_id, seq_id_target, resis, v=True):
        """
        :param seq_id: seq_id, 'AAML04000013.1/228868-228953'
        :param seq_id_target: seq_id of target, 'CP000721.1/2204691-2204778'
        :param resis: list resis, [5,6]

        map::

            [4, 5, 6]
            UAU-A
            UAU-AA
            UAU-AAC
            [5, 6, 7]
            CAC-U
            CAC-U-
            CAC-U-U
            [4, None, 5]

        """
        # get number for first seq
        print resis
        nresis = []
        for s in self.io:
            if s.id.strip() == seq_id.strip():
                    for i in resis:
                            print s.seq[:i+1] # UAU-A
                            nresis.append(i + s.seq[:i].count('-'))
        print nresis
        
        print self.map_seq_on_align(seq_id_target, nresis)
        return
    
        resis_target = []
        for s in self.io:
            if s.id.strip() == seq_id_target.strip():
                for i in nresis:
                    if v:print s.seq[:i]                    
                    if s.seq[i-1] == '-':
                        resis_target.append(None)
                    else:
                        resis_target.append(i - s.seq[:i].count('-'))
        return resis_target

    def map_seq_on_align(self, seq_id, resis,v=True):
        """
        :param seqid: seq_id, 'CP000721.1/2204691-2204775'
        :param resis: list resis, [5,6]

        maps::

            [5, 6, 8]
            CAC-U
            CAC-U-
            CAC-U-UA
            [4, None, 6]

        """
        if v: print resis
        nresis = []
        for s in self.io:
            if s.id.strip() == seq_id.strip():
                for i in resis:
                    if v:print s.seq[:i]                    
                    if s.seq[i-1] == '-':
                        nresis.append(None)
                    else:
                        nresis.append(i - s.seq[:i].count('-'))
        return nresis
    
    def head(self):
        return '\n'.join(self.lines[:5])

    def tail(self):
        return '\n'.join(self.lines[-5:])

    def describe(self):
        return str(self.io).split('\n')[0]
    
    def find_core(self, ids=[]):
        """Find common core for ids.

        :param id: list, ids of seq in the alignment to use
        """
        if ids == []:
            for s in self.io:
                ids.append(s.id)

        xx = list(range(0,len(self.io[0])))
        for i in range(0,len(self.io[0])): # if . don't use it
            for s in self.io:
                #print s.id
                if s.id in ids:
                    if s.seq[i] == '-':
                        xx[i] = '-'
                        break
                    else:
                        xx[i] = 'x'

        return ''.join(xx)
        shift = self.get_shift_seq_in_align()

        fnlist = open(self.fn).read().strip().split('\n')
        fnlist.insert(-2, 'x' + ' ' * (shift - 1) + ''.join(xx))
        #print fnlist
        for l in fnlist:
            print l

    def find_seq_in_align(self, seq, verbose=False):
        """
        :param seq: string, seq
        :param verbose: boolean, be verbose or not
        """
        seq = seq.replace('-','')
        for s in self.io:
            seq_str = str(s.seq).replace('-','')
            if verbose:
                print (seq_str)
            if seq == seq_str:
                print 'Match:', s.id
                print s
                print seq
                return s
        print 'Not found'

    def get_rfam_ss_notat_to_dot_bracket_notat(self, c):
        """Take (c)haracter and standardize ss"""
        if c in [',', '_']:
            return '.'
        if c == '<':
            return '('
        if c == '>':
            return ')'
        return c

    def get_clean_ss(self, ss):
        nss = ''
        for s in ss:
            nss += self.get_rfam_ss_notat_to_dot_bracket_notat(s)
        return nss
    
    def get_seq_ss(self,seq_id):#seq,ss):
        s = self.get_seq(seq_id).seq
        #print seq, ss
        # new
        nseq = ''
        nss = ''
        for i,j in zip(seq,ss):
            if i != '-': # gap
                #print i,j
                nseq += i
                nss += rfam_ss_notat_to_dot_bracket_notat(j)
        return nseq.strip(), nss.strip()

    def get_seq(self, seq_id):
        for s in self.io:
            if s.id == seq_id:
                return s
        raise Exception('Seq not found')

    def __str__(self):
        return (str(self.io))
    
if __name__ == '__main__':
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    #print a.get_shift_seq_in_align()
    #print a.map_seq_on_align('CP000721.1/2204691-2204778', [5,6,8])
    print a.map_seq_on_seq('AAML04000013.1/228868-228953', 'CP000721.1/2204691-2204778', [4,5,6])
    
    #print(record.letter_annotations['secondary_structure'])
    seq = a.get_seq('AAML04000013.1/228868-228953')
    print seq.seq
    print seq.ss
    print a.ss_cons

    print 'x'
    for s in a.io:
        print s.seq
        print s.ss_clean #letter_annotations['secondary_structure']
        
    for s in a.io:
        print s.seq_nogaps
        print s.ss_nogaps
