#./triplexibility2.py -t data/t1-1-Triple_cWW_cHS_GCU_exemplar_rpr_fx.pdb --result data/db_trx2.csv --sort --save
# python splix2.py ''
#apython trx_score_growth.py
#./triplexibility2.py -t data/gt11_cgC_rpr.pdb --result data/db_trx2.csv --sort --save

#python trx.py --edge cWW_cHS gcc
python trx_analysis.py 2>&1 | tee trx_analysis.txt
