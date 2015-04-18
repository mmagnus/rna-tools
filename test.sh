./yapdb_parser.py -h | tee yapdb_parser.out

./yapdb_parser.py -r input/1xjr.pdb 
./yapdb_parser.py -c input/1xjr.pdb > output/1xjr_clx.pdb

./yapdb_parser.py --rosetta2generic input/farna.pdb > output/farna_clx.pdb

./yapdb_parser.py --getchain A input/1xjr.pdb > output/1xjr_A_clx.pdb

./yapdb_parser.py --getseq input/1xjr.pdb > output/1xjr.seq
