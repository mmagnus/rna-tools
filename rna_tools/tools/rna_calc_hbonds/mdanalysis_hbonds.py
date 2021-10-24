import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

pdb = 'Triple_cWW_cHW_UAU_rpr_addChimera_onlyBases.pdb'
u = MDAnalysis.Universe(pdb)#s, pdb)

hbonds = HBA(universe=u)
hbonds.run()
