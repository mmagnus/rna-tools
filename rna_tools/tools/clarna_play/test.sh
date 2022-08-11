rm test_output/3e5fA_M15_c.pdb.outCR
rm test_output/3e5fA_M500_c.pdb.outCR

rna_clarna_run.py -ipdb test_data/3e5fA_M15_c.pdb | tee test_output/3e5fA_M15_c.pdb.outCR
rna_clarna_run.py -ipdb test_data/3e5fA_M500_c.pdb | tee test_output/3e5fA_M500_c.pdb.outCR
rna_clarna_run.py -ipdb test_data/3e5fA_M500_c.pdb | tee test_data/3e5fa_no_bps.pdb.outCR

rna_clarna_compare.py -iref test_output/3e5fA_M15_c.pdb.outCR -ichk test_output/3e5fA_M500_c.pdb.outCR
#lib/ClaRNAwd_to_vienaSS/ClaRNAwd_output_parser_get_SS test_output/3e5fA_M15_c.pdb.outCR
