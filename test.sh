./clarna_run.py -ipdb test_data/3e5fA_M15_c.pdb > test_output/3e5fA_M15_c.pdb.outCR
./clarna_run.py -ipdb test_data/3e5fA_M500_c.pdb > test_output/3e5fA_M500_c.pdb.outCR
./clarna_compare.py -iref test_output/3e5fA_M15_c.pdb.outCR -ichk test_output/3e5fA_M500_c.pdb.outCR
