echo rnastruc_clanstix.py
echo ------------------------------------------------------------
echo input test_data/matrix.txt
python rna_clanstix.py test_data/matrix.txt clans_run.txt
echo clans_run.txt generated

python rna_clanstix.py --groups-auto 10 --color-by-homolog --shape-by-source test_data/gmp_ref_mapping_pk_refX.txt test_data/gmp_color_by_homolog_shape_by_source.clans
python rna_clanstix.py --groups-auto 10 --color-by-homolog test_data/gmp_ref_mapping_pk_refX.txt test_data/gmp_color_by_homolog.clans
python rna_clanstix.py test_data/gmp_ref_mapping_pk_refX.txt test_data/gmp.clans
