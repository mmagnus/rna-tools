#!/usr/bin/env python
import os
from rna_pdb_tools.utils.rna_alignment.rna_alignment import RNAalignment

path = os.path.dirname(os.path.realpath(__file__))

alignment = RNAalignment(path + os.sep + 'test_data/RF00167.stockholm.sto')
print(alignment.tail())
print(alignment.ss_cons)
print(alignment[1])

alignment.plot()

from rna_pdb_tools.Seq import RNASequence
