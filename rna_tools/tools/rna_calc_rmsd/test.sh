set -x

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --ignore-files CGA --way backbone+sugar > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --save --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb  test_data/triples/*.pdb --save --suffix 'alignedXXX' --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

# ./rna_calc_rmsd_biopython.py -t test_data/triples/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/Triple_cWW_tSH_GCA_exemplar_rpr_ren.pdb  --way backbone+sugar --save --column-name 'AAAA' --triple-mode > test_output/2nd_triplex_FB_1AUA3_rpr_Triple_cWW_tSH_GCA_exemplar_rpr_ren_way_backbone+sugar_save_column-name_AAAA_triple-mode.txt

# ./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way c1p -o test_data/ways/rmsd_c1p.csv

# ./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way bases -o test_data/ways/rmsd_bases.csv

# ./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way all -o test_data/ways/rmsd_all.csv

# ./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr -sr -o test_data/ways/rmsd_all.csv

# ./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way backbone+sugar -o test_data/ways/backbone+sugar.csv

./rna_calc_rmsd.py -t test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb --target-selection A:1-47+52-62 --model-selection A:1-47+52-62 --model-ignore-selection A/57/O2\'+A/58/O2\' test_data/pistol/clusters/*.pdb test_data/pistol/clusters/pistol_thrs0.50A_clust01-000001_AA.pdb --target-column-name -o test_output/pistol_rmsd.csv

./rna_calc_rmsd.py -t test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb --target-selection A:1-47+52-62 --model-selection A:1-47+52-62 --model-ignore-selection A/57/O2\'+A/58/O2\' test_data/pistol/clusters/*.pdb test_data/pistol/clusters/pistol_thrs0.50A_clust01-000001_AA.pdb --name-rmsd-column RMSDtoPistol -o test_output/pistol_rmsd_name-rmsd-column.csv

# error ;-) fixed with the next line
##./rna_calc_rmsd_multi_targets.py --models test_data/multi-targets/rp21/*.pdb --targets test_data/multi-targets/rp21/solutions/*.pdb
#./rna_calc_rmsd_multi_targets.py --models test_data/multi-targets/rp21/*.pdb --targets test_data/multi-targets/rp21/solutions/*.pdb --target-selection A:1-27+29-41 --model-selection A:1-27+29-41
