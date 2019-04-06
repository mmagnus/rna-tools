#!/bin/bash

python renum_pdb_to_aln.py --residue_index_start 1 obj1 test_data/ALN_OBJ1_OBJ2.fa test_data/obj01.pdb
head test_data/obj01.pdb
echo
head test_data/obj01_out.pdb

tail test_data/obj01.pdb
echo
tail test_data/obj01_out.pdb


python renum_pdb_to_aln.py --residue_index_start 46 obj1 test_data/ALN_OBJ1_OBJ2.fa test_data/obj01.pdb 
head test_data/obj01.pdb
echo
head test_data/obj01_out.pdb

tail test_data/obj01.pdb
echo
tail test_data/obj01_out.pdb
