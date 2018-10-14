#!/usr/bin/python
#-*- coding: utf-8 -*-
from rna_pdb_tools.Seq import RNASequence
from rna_pdb_tools.rpt_config import CONTEXTFOLD_PATH, RNASTRUCTURE_PATH

seq = RNASequence("CCuuuuGG")
print(seq.predict_ss("mcfold", verbose=True))
