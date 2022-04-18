#./rna_rosetta_cluster.py test_data/selected.out 10
#./rna_rosetta_run.py -r -g -n 10 -c 1 --sandbox test_data test_data/tetraloop.fa
#python rna_rosetta_silent_split.py -o test_data test_data/selected.out
./rna_rosetta_get_score.py test_data/selected.out
./rna_rosetta_silent_random.py 5 test_data/selected.out
./rna_rosetta_silent_random.py --keep-order 5 test_data/selected.out
diff test_data/selected.out test_data/selected.out.rand5.out # should be nothing

