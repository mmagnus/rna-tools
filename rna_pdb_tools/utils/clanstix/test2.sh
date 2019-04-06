echo rnastruc_clanstix.py
echo ------------------------------------------------------------
echo input test_data/matrix.txt
python rna_clanstix.py --groups 10+10+1 test_data/matrix.txt clans_run.txt
echo clans_run.txt generated
cat clans_run.txt
