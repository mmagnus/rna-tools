#!/bin/bash

set -x
./rna_refinement.py test_data/input.pdb -s 10
./rna_refinement.py test_data/input.pdb -o test_data/output.pdb -s 10
./rna_refinement.py test_data/rp17_179c48aa-c0d3-4bd6-8e06-12081da22998_ALL_thrs6.20A_clust01-000001_AA._refx.pdb -s 10
rm -r tmp
