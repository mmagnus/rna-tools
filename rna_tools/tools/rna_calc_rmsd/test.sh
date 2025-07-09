#pymol -cq rna_calc_rmsd.py -- args /Users/magnus/Desktop/rmsd-all-vs-all/6qn3_A.cif  /Users/magnus/Desktop/rmsd-all-vs-all/5mqf_Z.cif

#./rna_calc_rmsd.py -m align -t /Users/magnus/Desktop/rmsd-all-vs-all/6qn3_A.cif  /Users/magnus/Desktop/rmsd-all-vs-all/5mqf_Z.cif
/Applications/PyMOL3.app/Contents/bin/python3 rna_calc_rmsd_pymol.py -m align -t /Users/magnus/Desktop/rmsd-all-vs-all/6qn3_A.cif  /Users/magnus/Desktop/rmsd-all-vs-all/5mqf_Z.cif
cat rmsds.csv
#
#set -x
#
#p=../../input/comparison
#./rna_calc_rmsd.py -t $p/4GXY_min.pdb $p/4GXY_min_reconstruction.pdb
