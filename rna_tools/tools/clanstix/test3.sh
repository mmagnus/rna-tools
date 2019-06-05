python rna_clanstix.py --multiprocessing test_data/gmp_ref_mapping_pk_refX.txt \
       --output test_data/gmp_mp.clans \
       --output-pmatrix \
       --output-pmatrix-fn test_data/foo_mp.txt

python rna_clanstix.py test_data/gmp_ref_mapping_pk_refX.txt \
       --output-pmatrix --output test_data/gmp_mp.clans \
       --output-pmatrix-fn test_data/foo.txt
