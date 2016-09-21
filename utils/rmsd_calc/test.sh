cmd="./rmsd_calc_biopy_dir.py -i test_data -o test_output/rmsd_calc_biopy_dir_matrix.tsv"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_dir.py -i test_data -o test_output/rmsd_calc_dir.tsv"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_to_target.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_to_target.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*"
echo $cmd
$cmd
echo

cmd="./rmsd_calc_to_target.py -t test_data/struc1.pdb --target_selection A:10 --model_selection A:10 --model_ignore_selection A/10/O2\' --target_ignore_selection A/10/O2\' -o test_output/rmsd_calc_dir_to_target.tsv test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb"
echo $cmd
$cmd
