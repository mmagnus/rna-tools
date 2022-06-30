./rna_df_merge.py test_data/*rmsd.csv --drop-col fn -o test_data/rna_df_merge.csv
./rna_df_merge_on.py fn test_data/min_*  -o test_data/rna_df_merge_on.csv
./rna_df_concat.py test_data/d3.csv test_data/d8.csv -o test_data/rna_df_concat.py
