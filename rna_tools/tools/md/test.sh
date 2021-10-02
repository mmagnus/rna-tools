set -x
# test for 2 chains
rna_gromacs_ready.py test_data/RA5_G.pdb | tee test_output/RA5_G.pdb
rna_minimize.py test_data/1duq_min_rpr_reduce_chain1.pdb | tee test_output/1duq_min_rpr_reduce_chain1.log
rna_get_energy.py test_data/1duq_min_rpr_reduce_chain1.pdb | tee test_output/1duq_min_rpr_reduce_chain1.log
