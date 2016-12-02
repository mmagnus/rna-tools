#./rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -o test_output/pistol_inf.csv -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb

#((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))'
./rna_calc_inf.py -s '((((([[[[[[)))))........(.((....(]]]]]].)..(((......)))...)).)'  -o test_output/pistol_inf.csv -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb ../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb
