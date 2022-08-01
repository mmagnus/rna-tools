rm test_output/*
./rna_mq_farfar2.py ../test/1y26X_output4_01-000001_AA.pdb | tee test_output/ff2.csv # && open test/farna.csv
./rna_mq_farfar2.py -r ../test/1y26X_output4_01-000001_AA.pdb | tee test_output/ff2_hires.csv # && open test/farna_hires.csv

../rna_mq_collect.py -f ../test/1y26X_output4_01-000001_AA.pdb -t FARFAR2_hires -o test_output/ff2_hires_rna_mq_collect.csv
../rna_mq_collect.py -f ../test/1y26X_output4_01-000001_AA.pdb -t FARFAR2 -o test_output/ff2_rna_mq_collect.csv
