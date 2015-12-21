./pdb_parser_lib.py

./yapdb_parser.py -h | tee yapdb_parser.out

./yapdb_parser.py -r input/1xjr.pdb 

./yapdb_parser.py --nohr -c input/1xjr.pdb > output/1xjr_clx.pdb

./yapdb_parser.py --nohr --rosetta2generic input/farna.pdb > output/farna_clx.pdb

./yapdb_parser.py --getchain A input/1xjr.pdb > output/1xjr_A_clx.pdb

./yapdb_parser.py --getseq input/1xjr.pdb > output/1xjr.seq

./yapdb_parser.py --nohr --clean input/1a9l_NMR_1_2_models.pdb > output/1a9l_NMR_1_2_models_tool.pdb

./yapdb_parser.py --nohr --clean input/1xjr_GTP.pdb > output/1xjr_GTP.pdb

./yapdb_parser.py --nohr --getrnapuzzle input/1xjr_onlyGTP.pdb

./yapdb_parser.py --nohr --getrnapuzzle input/1xjr_onlyGTP.pdb > input/1xjr_onlyGTP_rnapuzzle_ready.pdb

./yapdb_parser.py --nohr --clean input/1osw.pdb > output/1osw_NMR_1.pdb

./yapdb_parser.py --nohr --getrnapuzzle input/4GXY_3firstNt.pdb

./yapdb_parser.py --nohr --getrnapuzzle input/gtp.pdb  > output/gtp.pdb

./yapdb_parser.py --nohr --getrnapuzzle input/377D.pdb # should finish with error

./yapdb_parser.py --nohr --get_simrna_ready input/1xjr_no_op3.pdb > output/1xjr_no_op3_simrna_ready.pdb
