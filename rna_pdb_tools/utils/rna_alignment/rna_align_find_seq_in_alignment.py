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

In input file and input alignment all `-` are removed.

Or use -v::

    $ rna_align_find_seq_in_alignment.py -f 4lvv.seq -a RF01831.stockholm.sto -v
    UCAGAGUAGAAAACGAUGCGUUAAGUGUCCAGCAGACGGGGAGUUGCUGCCGGAACGAAAAGCAAAGCUUGCGGUAUCGUUUUCGCAUCCCGCUGA
    GCAGAGUAGACACAUGUGCGUUAAGUGCCGGAUGAACAGGGAGUUGUCACCCGGACGAAAAGAAAAUCUUGCGGUACAUGAGUCGCAUCCCGCUGC
    GCAGAGUAGGUUUGUGUGCGUUAAGUGCUGGUUGAACAGGGAGUUGUCAGCCGGACGAAAAGAUUUUCUUGCGGUACACGAAUCGCAUCCCGCUGC

so you can combine this with `grep`::

   $ rna_align_find_seq_in_alignment.py -f 4lvv.seq -a RF01831.stockholm.sto -v | grep UGUUGGGGUAGGAAUUACCGAGUUAUUGUCCAGCGGAACGGAAUGUGAACGCU
   UGUUGGGGUAGGAAUUACCGAGUUAUUGUCCAGCGGAACGGAAUGUGAACGCUGGAAGGAAGUAUUAGCUGCUUGAUAAUUUCGCAUUCACCACA

.. warning:: requires http://biopython.org/DIST/docs/api/Bio.AlignIO.StockholmIO-module.html

"""

from rna_alignment import RNAalignment

import sys
import argparse
from Bio import AlignIO

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-a', '--alignment', help="alignment in stockholm format", required=True)
    parser.add_argument('-f', '--file', help="fasta seq", required=True)
    parser.add_argument('-v', '--verbose', help="", action='store_true', default=False)
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    #args = parser.parse_args(['-f' , 'test_data/4lvv.seq', '-a', 'test_data/RF01831.stockholm.sto'])

    seq = open(args.file).readline().strip()

    a = RNAalignment(args.alignment)
    a.find_seq_in_align(seq, args.verbose)
