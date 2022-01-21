set -x
#python splix.py #  --d2 UA --growth 1 # wt U57:A27.A53
#python seqs.py AAG-0 # U57:A27.A53
#python seqs.py AUG-1 # suppresion
#python seqs.py ACG-0.5 # half-suppression
#python seqs.py AGG-0 #
#open splix.csv

# triplexibility.py
#./triplexibility.py -t data/still-wrong/t2-4.pdb data/still-wrong/0.90-Triple_tHW_tHH_ACA_rpr_renm.pdb
#./triplexibility.py --triple-mode --tseq UUU --way backbone+sugar -t data/still-wrong/t2-4.pdb data/still-wrong/* --save --debug

./triplexibility.py --triple-mode --tseq aca --way backbone+sugar -t data/still-wrong/t2-4.pdb data/still-wrong/0.90-Triple_tHW_tHH_ACA_rpr_renm_fx.pdb --save --debug
# still-wrong/*

#--verbose
#0.90-Triple_tHW_tHH_ACA_rpr_renm.pdb --save --triple-mode 
#./triplexibility.py --triple-mode --way c1+Nx -t data/still-wrong/t2-4.pdb data/still-wrong/*
#0.90-Triple_tHW_tHH_ACA_rpr_renm.pdb
