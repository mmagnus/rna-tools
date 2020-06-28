set -x
#((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))'
rna_calc_inf.py -s test_data/pistol.fa -o test_data/pistol_inf.csv -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb

rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb --no-stacking -pr -o test_data/pistol_inf_nostacking.csv > test_data/pistol_inf_nostacking.log

rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb -pr -o test_data/pistol_inf_nostacking.csv > test_data/pistol_inf_nostacking.log

rna_calc_inf.py --target-selection A:1-3+16-18 \
                --model-selection A:1-3+16-18 \
                -t ../rna_calc_rmsd/test_data/pistol/clusters/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb \
                -f ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb \
                -pr -o test_data/pistol_inf_nostacking.csv \
                --renumber-residues
