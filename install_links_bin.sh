#!/bin/bash

rm bin/*
curr_dir=$(pwd)
ln -s $curr_dir/rna_pdb_tools/utils/rmsd_calc/rmsd_calc_to_target.py $curr_dir/bin/rmsd_calc_to_target.py
ln -s $curr_dir/rna_pdb_tools/rna-pdb-tools.py $curr_dir/bin/rna-pdb-tools.py
ln -s $curr_dir/rna_pdb_tools/utils/diffpdb/diffpdb.py $curr_dir/bin/diffpdb
ln -s $curr_dir/rna_pdb_tools/utils/rna_multimodels/rna-pdb-merge-into-one.py $curr_dir/bin/rna-pdb-merge-into-one.py
echo 'Installed in ./bin'
ls -l bin
