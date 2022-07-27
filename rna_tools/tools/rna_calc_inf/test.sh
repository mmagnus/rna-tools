set -x
#((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))'
rm ../rna_calc_rmsd/test_data/pistol/clusters/*sel.pdb
rna_calc_inf.py -s test_data/pistol.fa -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb -pr -o test_data/pistol_inf_from_SeqSS.csv 

rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb -m 10 -o test_data/pistol_inf_m10.csv -pr

rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb -m 1 -o test_data/pistol_inf_m1.csv -pr

# --no-stacking 
#rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb --no-stacking -pr -o test_data/pistol_inf_nostacking.csv -f

rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb -pr -o test_data/pistol_inf_stacking.csv -f

rna_calc_inf.py --target-selection A:1-3+16-18 \
                --model-selection A:1-3+16-18 \
                -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb \
                -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb \
                -pr -o test_data/pistol_inf_selection.csv \
                --renumber-residues
