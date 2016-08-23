#!/bin/bash

./pdb_parser_lib.py

./rna-pdb-tools.py -h | tee rna-pdb-tools.out

./rna-pdb-tools.py -r input/1xjr.pdb 

./rna-pdb-tools.py --no_hr -c input/1xjr.pdb > output/1xjr_clx.pdb

./rna-pdb-tools.py --no_hr --rosetta2generic input/farna.pdb > output/farna_clx.pdb

./rna-pdb-tools.py --get_chain A input/1xjr.pdb > output/1xjr_A_clx.pdb

./rna-pdb-tools.py --get_seq input/1xjr.pdb > output/1xjr.seq

./rna-pdb-tools.py --no_hr --clean input/1a9l_NMR_1_2_models.pdb > output/1a9l_NMR_1_2_models_tool.pdb

./rna-pdb-tools.py --no_hr --clean input/1xjr_GTP.pdb > output/1xjr_GTP.pdb

./rna-pdb-tools.py --no_hr --get_rnapuzzle input/1xjr_onlyGTP.pdb

./rna-pdb-tools.py --no_hr --get_rnapuzzle input/1xjr_onlyGTP.pdb > output/1xjr_onlyGTP_rnapuzzle_ready.pdb

./rna-pdb-tools.py --no_hr --clean input/1osw.pdb > output/1osw_NMR_1.pdb

./rna-pdb-tools.py --no_hr --get_rnapuzzle input/4GXY_3firstNt.pdb

./rna-pdb-tools.py --no_hr --get_rnapuzzle input/gtp.pdb  > output/gtp.pdb

./rna-pdb-tools.py --no_hr --get_rnapuzzle input/377D.pdb # should finish with error

./rna-pdb-tools.py --no_hr --get_simrna_ready input/1xjr_no_op3.pdb > output/1xjr_no_op3_simrna_ready.pdb

./rna-pdb-tools.py --delete A:10-60 input/rp17.out.1.pdb > output/rp17_rmA10-60.pdb

# ClashCalc
cd ./utils/ClashCalc/
./ClashCalc.py
