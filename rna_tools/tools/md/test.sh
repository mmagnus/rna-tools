set -x
# test for 2 chains
#rna_gromacs_ready.py test_data/RA5_G.pdb | tee test_output/RA5_G.pdb
#rna_minimize.py test_data/1duq_min_rpr_reduce_chain1.pdb | tee test_output/1duq_min_rpr_reduce_chain1.log
#rna_get_energy.py test_data/1duq_min_rpr_reduce_chain1.pdb | tee test_output/1duq_min_rpr_reduce_chain1.log

python rna_traj_range.py test_data/AF-Q5TCX8-F1-model_v1_core_Ctrim_mdr_i316m_px_minpadd0.5_MD.pdb -s 1 -e 3
python rna_traj_range.py --rm-htm test_data/AF-Q5TCX8-F1-model_v1_core_Ctrim_mdr_i316m_px_minpadd0.5_MD.pdb -s 1 -e 3
python rna_traj_range.py --rm-htm test_data/AF-Q5TCX8-F1-model_v1_core_Ctrim_mdr_i316m_px_minpadd0.5_MD.pdb -s 3 -e 5
