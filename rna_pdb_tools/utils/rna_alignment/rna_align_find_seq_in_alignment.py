#!/usr/bin/env python

"""Usage::

    $ python rna_align_find_seq_in_alignment.py -a test_data/RF00167.stockholm.sto -f test_data/xx.fa
    Match: CP000903.1/3101736-3101822
    ID: CP000903.1/3101736-3101822
    Name: CP000903.1
    Description: CP000903.1/3101736-3101822
    Number of features: 0
    /start=3101736
    /end=3101822
    /accession=CP000903.1
    Seq('CAC-U-CGUAUAUACUCGGUAAUAUGG-UCCGAGC-GUUUCUACCUAGUUCCCA...UU-', SingleLetterAlphabet())

.. warning:: requires http://biopython.org/DIST/docs/api/Bio.AlignIO.StockholmIO-module.html

"""

import sys
import argparse
from Bio import AlignIO

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-a', '--alignment', help="alignment in stockholm format", required=True)
    parser.add_argument('-f', '--file', help="fasta seq", required=True)
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args() 

    seq = open(args.file).readline().strip()
    seq = seq.replace('-','')
    
    a = AlignIO.read(args.alignment, "stockholm")
    for s in a:
        seq_str = str(s.seq).replace('-','')
        if seq == seq_str:
            print 'Match:', s.id
            print s
            print seq
            sys.exit(0)
    print 'Not found'


