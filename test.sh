./pdb_parser_lib.py

./yapdb_parser.py -h | tee yapdb_parser.out

./yapdb_parser.py -r input/1xjr.pdb 

./yapdb_parser.py -c input/1xjr.pdb > output/1xjr_clx.pdb

./yapdb_parser.py --rosetta2generic input/farna.pdb > output/farna_clx.pdb

./yapdb_parser.py --getchain A input/1xjr.pdb > output/1xjr_A_clx.pdb

./yapdb_parser.py --getseq input/1xjr.pdb > output/1xjr.seq

./yapdb_parser.py -c input/1a9l_NMR_1_2_models.pdb > output/1a9l_NMR_1_2_models_tool.pdb

./yapdb_parser.py -c input/1xjr_GTP.pdb > output/1xjr_GTP.pdb

./yapdb_parser.py --getrnapuzzle input/1xjr_onlyGTP.pdb

./yapdb_parser.py --getrnapuzzle input/1xjr_onlyGTP.pdb > input/1xjr_onlyGTP_rnapuzzle_ready.pdb

./yapdb_parser.py --clean bucket/1osw.pdb > output/1osw_NMR_1.pdb
