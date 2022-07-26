#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
import fileinput

for line in fileinput.input():
    seq = Seq(line.strip().replace('-', ''))  # remove gaps
    protein = seq.translate()
    print(protein)
    protein_txt = ''
    index_txt = ''
    for index, p in enumerate(str(protein)):
         protein_txt += '\\' + p + '/'
         index_txt += '' + str(index + 1).ljust(3)
    print(protein_txt)
    print(index_txt)
