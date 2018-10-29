#!/usr/bin/python
#-*- coding: utf-8 -*-
from rna_pdb_tools.Seq import RNASequence

seq = RNASequence("CCuuuuGG")
print(seq.predict_ss("mcfold", verbose=True))

seq = RNASequence("CCuuuuGG")
print(seq.predict_ss("mcfold", constraints='(......)', verbose=True))
