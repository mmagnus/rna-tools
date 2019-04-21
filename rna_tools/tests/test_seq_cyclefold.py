#!/usr/bin/python

from rna_tools.Seq import RNASequence
from rna_tools.rna_tools_config import CONTEXTFOLD_PATH, RNASTRUCTURE_PATH

seq = RNASequence("CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAAAUAUCAGgUGCAA")
print(seq.predict_ss("rnastructure_CycleFold", verbose=True))
