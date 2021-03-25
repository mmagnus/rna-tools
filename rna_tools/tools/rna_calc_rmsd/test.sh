set -x 
./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --ignore-files CGA --way backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --save --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb  test_data/triples/*.pdb --save --suffix 'alignedXXX' --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

mv -v test_data/triples/*aligned*pdb test_data/triples/aligned

./rna_calc_rmsd_biopython.py -t test_data/triples2/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples2/Triple_cWW_tSH_GCA_exemplar_rpr_ren.pdb  --way backbone+sugar --save --column-name 'AAAA' --triple-mode > test_output/2nd_triplex_FB_1AUA3_rpr_Triple_cWW_tSH_GCA_exemplar_rpr_ren_way_backbone+sugar_save_column-name_AAAA_triple-mode.txt

./rna_calc_rmsd_biopython.py -t test_data/triples2/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples2/*rpr.pdb  --way backbone+sugar --save --triple-mode > test_output/2nd_triplex_FB_1AUA3_rpr_Triple_cWW_tSH_GCA_exemplar_rpr_ren_way_backbone+sugar_save_column-name_triple-mode.txt

./rna_calc_rmsd_all_vs_all.py -i test_data -o test_output/rmsd_calc_dir.tsv
./python rna_calc_rmsd_all_vs_all.py -i test_data -o test_output/rmsd_calc_dir_align.mat -m align

./rna_calc_rmsd.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/*.pdb

./rna_calc_rmsd.py -m align -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target_align.tsv test_data/*.pdb

./rna_calc_rmsd.py -m fit -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target_fit.tsv test_data/*.pdb

./rna_calc_rmsd.py -t test_data/struc1.pdb --target-selection A:10 --model-selection A:10 --model-ignore_selection A/10/O2\' --target-ignore_selection A/10/O2\' -o test_output/rmsd_calc_dir_to_target.tsv test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb

./rna_calc_rmsd.py -t test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb --target-selection A:1-47+52-62 --model-selection A:1-47+52-62 --model-ignore-selection A/57/O2\' test_data/pistol/clusters/pistol_thrs0.50A_clust01-000001_AA.pdb test_data/pistol/clusters/*.pdb

./rna_calc_rmsd.py -t test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb --target-selection A:1-47+52-62 --model-selection A:1-47+52-62 --model-ignore-selection A/57/O2\'+A/58/O2\' test_data/pistol/clusters/*.pdb test_data/pistol/clusters/pistol_thrs0.50A_clust01-000001_AA.pdb 


./rna_calc_rmsd_multi_targets.py --models test_data/multi-targets/rp21/*.pdb --targets test_data/multi-targets/rp21/solutions/*.pdb

cd test_data
../rna_calc_rmsd_multi_targets.py --models multi-targets/rp21/*.pdb --targets multi-targets/rp21/solutions/*.pdb --target-selection A:1-27+29-41 --model-selection A:1-27+29-41

