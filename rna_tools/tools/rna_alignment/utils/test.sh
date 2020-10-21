#!/bin/bash
set -v
./mfa-get-seq-for-species.py test_data/RF00026_demo.fa Mustela test_data/mustela.fa
./mfa-get-seq-for-species.py --format test_data/RF00026_demo.fa Mustela test_data/mustela_format.fa
./rna_align_strip_stk.py test_data/mutant_list.stk > test_data/mutant_list.txt

./rna_alignment_get_species.py test_data/u5_rfam_u5only.stk
./rna_alignment_get_species.py test_data/u5_rfam_u5only.stk --one
