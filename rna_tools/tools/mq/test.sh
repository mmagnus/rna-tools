#rm test/dfire.csv
set -x
#rm test/rasp.csv  # RASP,Dfire
# -t AnalyzeGeometry,ClashScore
#'ClashScore', 'AnalyzeGeometry', 'SimRNA_0',   'RNAscore', 'eSCORE','RNAkb',
#'RASP', 'RNAkb_all', 'RNA3DCNN', 'Dfire', 'FARNA', 'FARNA_hires', 'FARFAR2', 'FARFAR2_hires'
#Dfire
# RNA3DCNN
rm test/mq.csv
"""
x Dfire
"""
m='eSCORE'
#./rna_mq_collect.py -f -t $m test/1xjrA*.pdb -o test/mq.csv -m 8 # -v #
./rna_mq_collect.py -f test/1xjrA*.pdb -o test/mq.csv -m 8 # -v #
cat test/mq.csv 
# ,AnalyzeGeometry
#./mqaprnascore.py -f -t 'ClashScore' test/1xjrA_M1.pdb.pdb -o test/all.csv -v
#./rna_mq_remote.py x
