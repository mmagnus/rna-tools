set -x
#rna_minimize.py Triple_cWW_cHW_UAU_rpr.pdb
#Triple_cWW_cHW_UAU_rpr_addChimera_onlyBases.pdb #
#Triple_cWW_cHW_UAU_rpr_addChimera.pdb # Triple_cWW_cHW_UAU_rpr_bases_only_addh.pdb
#Triple_cWW_cHW_UAU_rpr_addChimera_onlyBases.pdb

for f in test_data/7jnh-AUA_ha.pdb # Triple_cWW_cHW_UAU_rpr_addChimera_onlyBases.pdb
	 do
	     python pytraj_hbond.py $f
	     python pytraj_hbond.py $f
	     python mdtraj_hbonds.py $f
	     python rna_calc_hbonds.py $f
         done
