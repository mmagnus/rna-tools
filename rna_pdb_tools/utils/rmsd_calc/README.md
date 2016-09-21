rmsd_calc_dir.py
-------------------------------------------------------------------------------

The program calculates all-atom rmsds (root-mean-square deviation) for all PDB structures in a given folder and save them to a file, as a matrix (which can be read by <https://github.com/m4rx9/rnastruc_clanstix>).

Usage:

	$ ./rmsd_calc_dir.py -i test_data -o test_output/matrix.txt

Output:

    $ head test_output/matrix.txt
    #smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000001_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000002_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000003_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000004_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000005_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000006_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000007_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000008_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000009_AA.pdb smk_ready_samsoaked_chainB_fix.pdb-000001_chainA_allmove_smk_ready_samsoaked_chainA_SAM_fix.pdb-000001_allmove_cst_100_01_1-000010_AA.pdb 
	0.0 4.359 5.81 5.301 10.851 14.145 14.833 15.55 14.988 15.201 
	4.359 0.0 2.576 2.497 8.345 13.145 13.857 14.902 14.549 14.946 
	5.81 2.576 0.0 2.353 7.348 12.595 13.326 14.481 14.182 14.638 
	5.301 2.497 2.353 0.0 7.116 11.929 12.6 13.699 13.535 13.867 
	10.851 8.345 7.348 7.116 0.0 7.299 7.959 9.606 9.804 10.177 
	14.145 13.145 12.595 11.929 7.299 0.0 2.034 3.939 5.029 4.176 
	14.833 13.857 13.326 12.6 7.959 2.034 0.0 3.021 4.656 3.951 
	15.55 14.902 14.481 13.699 9.606 3.939 3.021 0.0 3.36 3.212 
	14.988 14.549 14.182 13.535 9.804 5.029 4.656 3.36 0.0 3.251 

The program is using Biopython:

Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878

rmsd_calc_to_target.py
-------------------------------------------------------------------------------

	./rmsd_calc_to_target.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb
	rmsd_calc_rmsd_to_target
	--------------------------------------------------------------------------------
	# of models: 4
	target:struc1.pdb      	rmsd_all
	struc1.pdb     	0.0
	struc2.pdb     	11.803
	struc3.pdb     	4.879
	struc4.pdb     	3.982
	tsv was created!  test_output/rmsd_calc_dir_to_target.tsv

Install
==========================
Add to your PATH in .bashrc something like this `/home/magnus/src/rna-pdb-tools/utils/rmsd_calc/`
