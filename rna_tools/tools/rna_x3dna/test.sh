set -x
./rna_x3dna.py test_data/1xjr.pdb > test_data/1xjr_output.pdb
./rna_x3dna.py --show-log test_data/CG.pdb > test_data/CG.log
