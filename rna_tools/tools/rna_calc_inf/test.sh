#./rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -o test_output/pistol_inf.csv -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb

#((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))'
./rna_calc_inf.py -s pistol.ss -o test_output/pistol_inf.csv -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb
#../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb

./rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb --no-stacking -pr -o test_output/pistol_inf_nostacking.csv > test_output/pistol_inf_nostacking.log

./rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb -pr -o test_output/pistol_inf_stacking.csv > test_output/pistol_inf_stacking.log
