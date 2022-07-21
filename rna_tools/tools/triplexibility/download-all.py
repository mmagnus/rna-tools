#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import os
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    txt = """AA	AAA AAC	AAG	AAU
AC	ACA	ACC	ACG	ACU
AG	AGA	AGC	AGG	AGU
AU	AUA	AUC	AUG	AUU
CA	CAA	CAC	CAG	CAU
CC	CCA	CCC	CCG	CCU
CG	CGA	CGC	CGG	CGU
CU	CUA	CUC	CUG	CUU
GA	GAA	GAC	GAG	GAU
GC	GCA	GCC	GCG	GCU
GG	GGA	GGC	GGG	GGU
GU	GUA	GUC	GUG	GUU
UA	UAA	UAC	UAG	UAU
UC	UCA	UCC	UCG	UCU
UG	UGA	UGC	UGG	UGU
UU	UUA	UUC	UUG	UUU
"""
    for l in txt.split('\n'):
        l = l.strip()
        triples = l.split()
        for t in triples[1:]:
            print(t)
            os.system('wget http://rna.bgsu.edu/triples/seq/' + t + '.html')
