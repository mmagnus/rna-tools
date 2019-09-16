#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import fileinput

for line in fileinput.input():
    seq = Seq(line.strip())
    protein = seq.translate()
    print(protein)
    protein_txt = ''
    for p in str(protein):
         protein_txt += '\\' + p + '/'
    print(protein_txt)
