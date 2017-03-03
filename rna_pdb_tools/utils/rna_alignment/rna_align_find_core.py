#!/usr/bin/env python

"""http://biopython.org/wiki/AlignIO
"""

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import rna_alignment as ra
reload(ra)

if __name__ == '__main__':
    fn = 'test_data/RF00167.stockholm.sto'
    ids = ['AL591975.1/251136-251218','CP000721.1/2204691-2204778']#ACCL02000010.1/116901-116991']
    
    a = ra.RNAalignment(fn)
    a.find_core(ids)
    
