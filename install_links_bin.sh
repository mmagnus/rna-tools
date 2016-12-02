#!/bin/bash

mkdir -p bin
rm bin/*
curr_dir=$(pwd)
ln -s $curr_dir/rna_pdb_tools/utils/rna_calc_rmsd/rna_calc_rmsd.py $curr_dir/bin/rna_calc_rmsd.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_calc_rmsd/rna_calc_rmsd_all_vs_all.py $curr_dir/bin/rna_calc_rmsd_all_vs_all.py
ln -s $curr_dir/rna_pdb_tools/rna-pdb-tools.py $curr_dir/bin/rna-pdb-tools.py
ln -s $curr_dir/rna_pdb_tools/utils/diffpdb/diffpdb.py $curr_dir/bin/diffpdb
ln -s $curr_dir/rna_pdb_tools/utils/rna_multimodels/rna-pdb-merge-into-one.py $curr_dir/bin/rna-pdb-merge-into-one.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_calc_inf/rna_calc_inf.py $curr_dir/bin/rna_calc_inf.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_convert_pseudoknot_formats/rna_pk_simrna_to_one_line.py $curr_dir/bin/rna_pk_simrna_to_one_line.py
ln -s $curr_dir/rna_pdb_tools/utils/clarna_app/clarna_app.py $curr_dir/bin/clarna_app.py
echo 'Installed in ./bin'
ls -l bin
echo
echo 'Broken links:'
find bin -type l -exec sh -c "file -b {} | grep -q ^broken" \; -print
echo '^ should be none!'
