#!/bin/bash

rm bin/*
curr_dir=$(pwd)
ln -s $curr_dir/rna_pdb_tools/utils/rmsd_calc/rmsd_calc_to_target.py $curr_dir/bin/rmsd_calc_to_target.py
echo 'Installed in ./bin'
ls bin
