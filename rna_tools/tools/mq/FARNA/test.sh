#parallel "rna_mq_farna.py {}" ::: test_data/*pdb | tee test_data/rna_mq_farna.csv
#parallel "rna_mq_farna.py -r {}" ::: test_data/*pdb | tee test_data/rna_mq_farna_hires_min.csv # minimize_cmd = ' -minimize_rna '
#rna_mq_farna.py -r test_data/*pdb | tee test_data/rna_mq_farna_hires_min.csv # minimize_cmd = ' -minimize_rna '
#rna_mq_farna.py -r test_data/1msy_output4_01-000001_AA_edited.pdb #| tee test_data/1_rna_mq_farna_hires_min.csv # minimize_cmd = ' -minimize_rna '
set -x
s=1y26X_output4_01-000001_AA.pdb
rm test_data/${s}*csv ############# !!!!!!!!!!!!!!!!!!!!! remove csv file !!!!!!!!!!!!!!!!!!!!
rna_mq_farna.py test_data/$s | tee test_data/${s}_rna_mq_collect_rna_mq_farna_lores.csv
rna_mq_collect.py -vf -t FARNA test_data/$s -o test_data/${s}_rna_mq_collecttFARNA_lores.csv # appends
rm test_data/*log
# /Users/magnus/Desktop/mqapRNAdb/decoys/curr/1y26A_simrna_unfolding/rna_csv.py test_data/${s}_rna_mq_collect_rna_mq_farna_lores.csv --flat
# hires
rna_mq_farna.py -r test_data/$s | tee test_data/${s}_rna_mq_collect_rna_mq_farna.csv
rna_mq_collect.py -vf -t FARNA_hires test_data/$s -o test_data/${s}_rna_mq_collecttFARNA_hires.csv

#cat test_data/1y26X_output4_01-000001_AA_rna_mq_collecttFARNA_hires*
#1msy_output4_01-000001_AA_edited.pdb #| tee test_data/1_rna_mq_farna_hires_min.csv # minimize_cmd = ' -minimize_rna '

