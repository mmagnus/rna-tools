#!/usr/bin/env python

"""Usage::

    $ rna_align_find_seq_in_alignment.py -a pistol_noPk.sto -s hcf.fa
    ('Match:', 'HCF12C_58327/229-301')
    ID: HCF12C_58327/229-301
    Name: HCF12C_58327
    Description: HCF12C_58327/229-301
    Number of features: 0
    /start=229
    /end=301
    /accession=HCF12C_58327
    Per letter annotation for: secondary_structure
    Seq('UGCCGUUUG-AGCGGCA-UUAAACA--GGUCU-UAAGCUCAA-AGCG-UCACCG...ACA', SingleLetterAlphabet())

In input file and input alignment all `-` are removed.

Or use -v::

    $ rna_align_find_seq_in_alignment.py -s 4lvv.seq -a RF01831.stockholm.sto -v
    UCAGAGUAGAAAACGAUGCGUUAAGUGUCCAGCAGACGGGGAGUUGCUGCCGGAACGAAAAGCAAAGCUUGCGGUAUCGUUUUCGCAUCCCGCUGA
    GCAGAGUAGACACAUGUGCGUUAAGUGCCGGAUGAACAGGGAGUUGUCACCCGGACGAAAAGAAAAUCUUGCGGUACAUGAGUCGCAUCCCGCUGC
    GCAGAGUAGGUUUGUGUGCGUUAAGUGCUGGUUGAACAGGGAGUUGUCAGCCGGACGAAAAGAUUUUCUUGCGGUACACGAAUCGCAUCCCGCUGC

so you can combine this with `grep`::

   $ rna_align_find_seq_in_alignment.py -s 4lvv.seq -a RF01831.stockholm.sto -v | grep UGUUGGGGUAGGAAUUACCGAGUUAUUGUCCAGCGGAACGGAAUGUGAACGCU
   UGUUGGGGUAGGAAUUACCGAGUUAUUGUCCAGCGGAACGGAAUGUGAACGCUGGAAGGAAGUAUUAGCUGCUUGAUAAUUUCGCAUUCACCACA

.. warning:: requires http://biopython.org/DIST/docs/api/Bio.AlignIO.StockholmIO-module.html

"""

from rna_tools.tools.rna_alignment.rna_alignment import RNAalignment

import sys
import argparse
from Bio import AlignIO


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-a', '--alignment', help="alignment in stockholm format", required=True)
    parser.add_argument('-s', '--seq', help="fasta seq", required=True)
    parser.add_argument('-v', '--verbose', help="", action='store_true', default=False)
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    #args = parser.parse_args(['-f' , 'test_data/4lvv.seq', '-a', 'test_data/RF01831.stockholm.sto'])

    f = open(args.seq)
    header = f.readline().strip()
    seq = f.readline().strip()

    a = RNAalignment(args.alignment)
    a.find_seq(seq, args.verbose)
