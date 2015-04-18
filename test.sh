./yapdb_parser.py -h | tee yapdb_parser.out

./yapdb_parser.py -r input/1xjr.pdb
./yapdb_parser.py -c input/1xjr.pdb

./pdb_checker.py
