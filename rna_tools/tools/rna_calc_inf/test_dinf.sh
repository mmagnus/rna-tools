./rna_calc_dinf.py  test_data/1Y26.pdb > test_data/1Y26_dinf.pdb.outCR
./rna_calc_dinf.py  test_data/1y26X_M451.pdb > test_data1y26X_M451_dinf.pdb.outCR

clarna_compare.py -v -iref test_output/1y26X_M451_dinf.pdb.outCR -ichk test_output/1Y26_dinf.pdb.outCR | tee  test_output/1Y26_dinf_comparison.txt
