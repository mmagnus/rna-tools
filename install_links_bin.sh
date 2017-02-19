#!/bin/bash

mkdir -p bin
rm bin/*
curr_dir=$(pwd)
ln -s $curr_dir/rna_pdb_tools/utils/rna_calc_rmsd/rna_calc_rmsd.py $curr_dir/bin/rna_calc_rmsd.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_calc_rmsd/rna_calc_rmsd_all_vs_all.py $curr_dir/bin/rna_calc_rmsd_all_vs_all.py
ln -s $curr_dir/rna_pdb_tools/rna_pdb_tools.py $curr_dir/bin/rna_pdb_tools.py
ln -s $curr_dir/rna_pdb_tools/utils/diffpdb/diffpdb.py $curr_dir/bin/diffpdb
ln -s $curr_dir/rna_pdb_tools/utils/rna_multimodels/rna_pdb_merge_into_one.py $curr_dir/bin/rna_pdb_merge_into_one.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_calc_inf/rna_calc_inf.py $curr_dir/bin/rna_calc_inf.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_convert_pseudoknot_formats/rna_pk_simrna_to_one_line.py $curr_dir/bin/rna_pk_simrna_to_one_line.py
ln -s $curr_dir/rna_pdb_tools/utils/clarna_app/clarna_app.py $curr_dir/bin/clarna_app.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_helix_vis/rna_helix_vis.py $curr_dir/bin/rna_helix_vis.py
ln -s $curr_dir/rna_pdb_tools/utils/misc/rna_add_chain.py $curr_dir/bin/rna_add_chain.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_sali2dotbracket/rna_sali2dotbracket.py $curr_dir/bin/rna_sali2dotbracket.py

ln -s $curr_dir/rna_pdb_tools/utils/rna_rosetta/rna_rosetta_run.py $curr_dir/bin/rna_rosetta_run.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_rosetta/rna_rosetta_cluster.py $curr_dir/bin/rna_rosetta_cluster.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_rosetta/rna_rosetta_min.py $curr_dir/bin/rna_rosetta_min.py
ln -s $curr_dir/rna_pdb_tools/utils/rna_rosetta/rna_rosetta_n.py $curr_dir/bin/rna_rosetta_n.py
ln -s $curr_dir/rna_pdb_tools/utils/simrna_trajectory/rna_simrna_lowest.py $curr_dir/bin/rna_simrna_lowest.py
ln -s $curr_dir/rna_pdb_tools/utils/simrna_trajectory/rna_simrna_extract.py $curr_dir/bin/rna_simrna_extract.py

echo 'Installed in ./bin'
ls -l bin
echo
echo 'Broken links'
echo '-- below this should be none --'
find bin -type l -exec sh -c "file -b {} | grep -q ^broken" \; -print
echo '-- ^ should be none! -- '
