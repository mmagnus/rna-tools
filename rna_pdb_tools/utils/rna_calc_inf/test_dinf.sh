./rna_calc_dinf.py  test_output/1Y26.pdb > test_output/1Y26_dinf
./rna_calc_dinf.py  test_output/1y26X_M451.pdb > test_output/1y26X_M451_dinf

clarna_compare.py -iref test_output/1y26X_M451_dinf.pdb.outCR -ichk test_output/1Y26_dinf.pdb.outCR | tee  test_output/1Y26_dinf_comparison.txt
