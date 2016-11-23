rna_calc_inf
-------------------------------------------------------------------------------

	$ rna_calc_inf.py -t ../rmsd_calc/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb \
		../rmsd_calc/test_data/pistol/clusters/*.pdb \
		-o test_output/pistol_inf.csv

	rna-calc-inf
	--------------------------------------------------------------------------------
	target, fn, inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC
	5k7c_clean_onechain_renumber_as_puzzle_srr.pdb.outCR     pistol_thrs0.50A_clust01-000001_AA.pdb.outCR      0.642      NA         0.874      0.000      0.944      0.810      0.000      0.000
	5k7c_clean_onechain_renumber_as_puzzle_srr.pdb.outCR     pistol_thrs0.50A_clust02-000001_AA.pdb.outCR      0.642      NA         0.865      0.000      0.889      0.842      0.000      0.000
	5k7c_clean_onechain_renumber_as_puzzle_srr.pdb.outCR     pistol_thrs0.50A_clust03-000001_AA.pdb.outCR      0.577      NA         0.754      0.000      0.833      0.682      0.000      0.000
	csv was created!  test_output/pistol_inf.csv

ClaRNA_play required!
https://gitlab.genesilico.pl/RNA/ClaRNA_play (internal GS gitlab server)
