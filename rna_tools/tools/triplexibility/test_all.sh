./triplexibility2.py -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --result data/db_trx2.csv --sort --save

#set -x
#python splix.py #  --d2 UA --growth 1 # wt U57:A27.A53
#python seqs.py AAG-0 # U57:A27.A53
#python seqs.py AUG-1 # suppresion
#python seqs.py ACG-0.5 # half-suppression
#python seqs.py AGG-0 #
#open splix.csv

# triplexibility.py
#./triplexibility.py -t data/still-wrong/t2-4.pdb data/still-wrong/0.90-Triple_tHW_tHH_ACA_rpr_renm.pdb
#./triplexibility.py --triple-mode --tseq UUU --way backbone+sugar -t data/still-wrong/t2-4.pdb data/still-wrong/* --save --debug


#./triplexibility.py --triple-mode --tseq aca --way backbone+sugar -t data/still-wrong/t2-4.pdb data/still-wrong/0.90-Triple_tHW_tHH_ACA_rpr_renm_fx.pdb --save --debug
# still-wrong/*

#--verbose
#0.90-Triple_tHW_tHH_ACA_rpr_renm.pdb --save --triple-mode 
#./triplexibility.py --triple-mode --way c1+Nx -t data/still-wrong/t2-4.pdb data/still-wrong/*
#0.90-Triple_tHW_tHH_ACA_rpr_renm.pdb

# --tseq aca #  --way backbone+sugar
# data/still-wrong/t2-4.pdb 
# --way backbone+sugar#
#./triplexibility.py --triple-mode -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb db/triples-all-v2-rpr/* --result data/spl-1-1.csv --sort #--save #--debug
./triplexibility.py --triple-mode -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --files db/triples-all-v2-rpr/* --result data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.csv --sort #--save #--debug

./triplexibility.py --triple-mode -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --files \
		   db/triples-all-v2-rpr/*cWW_cHS_GCU* --result data/cWW_cHS_GCU.csv --sort #--save #--debug

./triplexibility.py --triple-mode -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --files \
		   db/triples-all-v2-rpr/*cWW* --result data/XcWWX.csv --sort #--save #--debug

./triplexibility.py --triple-mode -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --files \
		   db/triples-all-v2-rpr/*cWW* --result data/XcWWX.csv --sort --save #--debug

./triplexibility.py --triple-mode -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --files \
		   db/triples-all-v2-rpr/*cWW* --result data/XcWWX.csv --sort --save #--debug

./triplexibility.py --triple-mode -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --result data/db.csv --sort # --save #--debug

