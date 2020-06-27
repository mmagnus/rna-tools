rna_calc_inf
-------------------------------------------------------------------------------

````
rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb \
                ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb \
                -o test_data/pistol_inf.csv --print-results
100% (4 of 4) |#########################################################################################################################################| Elapsed Time: 0:00:00 ETA:  00:00:00
csv was created!  test_data/pistol_inf.csv
                                           target                                              fn  inf_all  inf_stack  inf_WC  inf_nWC  sns_WC  ppv_WC  sns_nWC  ppv_nWC
0  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb          pistol_thrs0.50A_clust02-000001_AA.pdb     0.41        0.0    0.86     0.38    0.89    0.84     0.33     0.43
1  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb          pistol_thrs0.50A_clust01-000001_AA.pdb     0.42        0.0    0.87     0.35    0.94    0.81     0.33     0.38
2  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb     0.57        0.0    1.00     1.00    1.00    1.00     1.00     1.00
3  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb          pistol_thrs0.50A_clust03-000001_AA.pdb     0.37        0.0    0.75     0.30    0.83    0.68     0.22     0.40
```
