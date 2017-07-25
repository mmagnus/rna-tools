set -x
#--bfactor#
./rna_pdb_edit_occupancy_bfactor.py --occupancy --select A:1-2 --set-not-selected-to 20 --set-to 0.5 -o test_data/3w3s_homologymodel_out.pdb test_data/3w3s_homologymodel.pdb 
head -n 50 test_data/3w3s_homologymodel_out.pdb
