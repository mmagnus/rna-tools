./yapdb_parser.py -h | tee yapdb_parser.out

./yapdb_parser.py -r input/1xjr.pdb bug
./yapdb_parser.py -c input/1xjr.pdb output/1xjr_clx.pdb

./yapdb_parser.py --rosetta2generic input/farna.pdb output/farna_clx.pdb
