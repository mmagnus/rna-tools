cmd="./rmsd_calc_biopy_dir.py -i test_data -o test_output/rmsd_calc_biopy_dir_matrix.tsv"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_dir.py -i test_data -o test_output/rmsd_calc_dir.tsv"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_to_target.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*.pdb"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_to_target.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*.pdb"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_to_target.py -t test_data/struc1.pdb --target_selection A:10 --model_selection A:10 --model_ignore_selection A/10/O2\' --target_ignore_selection A/10/O2\' -o test_output/rmsd_calc_dir_to_target.tsv test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb"
echo $cmd
$cmd

cmd="./rmsd_calc_to_target.py -t test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb --target_selection A:1-47+52-62 --model_selection A:1-47+52-62 /Users/magnus/work/src/rna-pdb-tools/rna_pdb_tools/utils/rmsd_calc/test_data/pistol/clusters/pistol_thrs0.50A_clust01-000001_AA.pdb --model_ignore_selection A/57/O2\'  test_data/pistol/clusters/*_AA.pdb"
echo $cmd
$cmd
