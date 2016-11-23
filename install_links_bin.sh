#!/bin/bash

rm bin/*
curr_dir=$(pwd)
ln -s $curr_dir/rna_pdb_tools/utils/rmsd_calc/rmsd_calc_to_target.py $curr_dir/bin/rmsd_calc_to_target.py
ln -s $curr_dir/rna_pdb_tools/rna-pdb-tools.py $curr_dir/bin/rna-pdb-tools.py
ln -s $curr_dir/rna_pdb_tools/utils/diffpdb/diffpdb.py $curr_dir/bin/diffpdb
ln -s $curr_dir/rna_pdb_tools/utils/rna_multimodels/rna-pdb-merge-into-one.py $curr_dir/bin/rna-pdb-merge-into-one.py
ln -s $curr_dir/rna_pdb_tools/utils/rna-calc-inf/rna-calc-inf.py $curr_dir/bin/rna-calc-inf.py
echo 'Installed in ./bin'
ls -l bin
echo
echo 'Broken links:'
find bin -type l -exec sh -c "file -b {} | grep -q ^broken" \; -print
echo '^ should be none!'
