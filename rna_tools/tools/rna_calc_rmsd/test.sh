set -x

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --ignore-files CGA --way backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/*.pdb --save --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

./rna_calc_rmsd_biopython.py -t test_data/2nd_triplex_FB_1AUA3_rpr.pdb  test_data/triples/*.pdb --save --suffix 'alignedXXX' --way=backbone+sugar \
			     > test_output/2nd_triplex_FB_CBA_ignore_files_CGA.csv

# ./rna_calc_rmsd_biopython.py -t test_data/triples/2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples/Triple_cWW_tSH_GCA_exemplar_rpr_ren.pdb  --way backbone+sugar --save --column-name 'AAAA' --triple-mode > test_output/2nd_triplex_FB_1AUA3_rpr_Triple_cWW_tSH_GCA_exemplar_rpr_ren_way_backbone+sugar_save_column-name_AAAA_triple-mode.txt

./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way c1p -o test_data/ways/rmsd_c1p.csv

./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way bases -o test_data/ways/rmsd_bases.csv

./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way all -o test_data/ways/rmsd_all.csv

./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr -sr -o test_data/ways/rmsd_all.csv

./rna_calc_rmsd.py --model-selection='A:52+53+59+60+61+80+B:21+22+23' --target-selection='A:52+53+59+60+61+80+B:21+22+23' -t test_data/ways/yC_5LJ3_U2U6_core_mdrFx_onlyTriplex_rpr.pdb test_data/ways/*fixChains.pdb -pr --way backbone+sugar -o test_data/ways/backbone+sugar.csv
