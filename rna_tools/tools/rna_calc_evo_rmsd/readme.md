
	./rna_calc_evo_rmsd.py -a test_data/aln_nox.aln -t test_data/4qk8_cl.pdb test_data/4qlm_cl.pdb -v

	 target name not provided; using basename: 4qk8_cl
	Alignment records:
	  4qk8_cl GU---UGCCGAAU---CCGAAAG-GU-A-CGGAGGAACCG---CUUUUUG----GGGUUAAUCUGC---AGUGA---AGCUG-----CAGUAGGGAUA-----CCUUCUG---UCCCGCACCCGACAGCUAACUCCGGAGGCAAUA-AA--GGAAGGA
	  4qlm_cl -A---UCGCUGAACG--------------CGGGGGACCCA--------------GGGGGCGAAUCU-CUUCCGAA--AGGAA-----GAGUAGGGUUA-----CUCCUUCG--ACCCGAGCCCGUCAGCUAACCUCGCAAGCGUCCGAA--GGAGAA-
	 Inferred selector from gapless columns.
	Selector source: gapless_columns
	Selector (x-line): -x---xxxxxxxx----------------xxxxxxxxxxx--------------xxxxxxxxxxxx---xxxxx---xxxxx-----xxxxxxxxxxx-----xxxxxxx---xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx-xx--xxxxxx-
	target test_data/4qk8_cl.pdb
	 Selected residues for 4qk8_cl: [2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119]
	 mapping file not provided; using PDB basenames as alignment IDs
	 # of rnastruc : 1
	 rnastruc: ['4qlm_cl:4qlm_cl']
	 WARNING: if any of your PDB file is missing, check mapping!
	 Selected residues for 4qlm_cl: [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 101, 102, 103, 104, 105, 106, 107, 108]
			target        model  rmsd group_name
	0  4qk8_cl.pdb  4qlm_cl.pdb  6.11


	rna_calc_evo_rmsd.py -a test_data/trna.sto -t test_data/1ehz_std.pdb test_data/6Y2L_2_std.pdb -v
	/Users/magnus/miniconda3/bin/rna_calc_evo_rmsd.py:4: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
	  __import__('pkg_resources').require('rna-tools==3.24.post0.dev3')
	 target name not provided; using basename: 1ehz_std
	Alignment records:
	  6Y2L_2_std GCCCGGAUAGCUCAGUcGGUAGAGCAGGGGAUUGAAAAUCCCCGUGuCCUUGGUUCGAUUCCGAGUCCGGGCAcca
	  1ehz_std GCGGAUUUAGCUCAGUuGGGAGAGCGCCAGACUGAAGAUCUGGAGGuCCUGUGUUCGAUCCACAGAAUUCGCAcca
	Selector source: reference_annotation
	Selector (x-line): GgagauaUAGCucAgU.GGUAgaGCgucgGaCUuaaAAuCcgaagg.cgcgGGUUCgAaUCCcgcuaucucCa...
	target test_data/1ehz_std.pdb
	 Selected residues for 1ehz_std: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73]
	 mapping file not provided; using PDB basenames as alignment IDs
	 # of rnastruc : 1
	 rnastruc: ['6Y2L_2_std:6Y2L_2_std']
	 WARNING: if any of your PDB file is missing, check mapping!
	 Selected residues for 6Y2L_2_std: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73]
			 target           model   rmsd group_name
	0  1ehz_std.pdb  6Y2L_2_std.pdb  1.652

Automated mode, no alignment needed. Sequences are taken from the PDB files and searched against Rfam (`cmscan`); if the target and a model hit the same family, they are aligned to the family covariance model (`cmalign`) and the residues in the consensus (`#=GC RF`) columns are used for the RMSD. Requires Infernal and Rfam.cm (`--rfam_db` or `RFAM_DB_PATH`, see RfamAlign.py). Alignments are saved to `<output>_rfam/<family>.sto`:

	rna_calc_evo_rmsd.py --rfam_db Rfam.cm -t test_data/1ehz_std.pdb test_data/6Y2L_2_std.pdb

Use `--rfam_evalue 0.01` instead of Rfam gathering thresholds (`--cut_ga`, default). Models that do not share a family with the target are skipped.

Without Rfam, `--auto pairwise` aligns each model to the target with a pairwise sequence alignment (all aligned residues are used):

	rna_calc_evo_rmsd.py --auto pairwise -t test_data/1ehz_std.pdb test_data/6Y2L_2_std.pdb

	         target           model  rmsd  n_residues family group_name
	0  1ehz_std.pdb  6Y2L_2_std.pdb   3.3          76

Residues are paired according to the alignment columns (the n-th residue of a sequence in the alignment is the n-th nucleotide in the PDB file), not by the residue numbers in the PDB files.
