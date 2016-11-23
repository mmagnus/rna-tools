cmd="./rna_calc_rmsd_all_vs_all.py -i test_data -o test_output/rmsd_calc_dir.tsv"
echo $cmd
$cmd
echo

cmd="./rna_calc_rmsd.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*.pdb"
echo $cmd
$cmd
echo

cmd="./rna_calc_rmsd.py -m align -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*.pdb"
echo $cmd
$cmd
echo

cmd="./rna_calc_rmsd.py -m fit -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*.pdb"
echo $cmd
$cmd
echo

cmd="./rna_calc_rmsd.py -t test_data/struc1.pdb --target_selection A:10 --model_selection A:10 --model_ignore_selection A/10/O2\' --target_ignore_selection A/10/O2\' -o test_output/rmsd_calc_dir_to_target.tsv test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb"
echo $cmd
$cmd

cmd="./rna_calc_rmsd.py -t test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb --target_selection A:1-47+52-62 --model_selection A:1-47+52-62 test_data/pistol/clusters/pistol_thrs0.50A_clust01-000001_AA.pdb --model_ignore_selection A/57/O2' test_data/pistol/clusters/*.pdb"
echo $cmd
$cmd
