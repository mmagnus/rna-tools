rna_calc_rmsd.py
-------------------------------------------------------------------------------

Using PyMOL align:

	./rna_calc_rmsd.py -m align -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb
	rmsd_calc_rmsd_to_target
	--------------------------------------------------------------------------------
	method: align
	# of models: 4
	struc1.pdb 0.0 1321
	struc2.pdb 11.8030576706 1321
	struc3.pdb 4.87882900238 1321
	struc4.pdb 3.98158359528 1321
	# of atoms used: 1321
	csv was created!  test_output/rmsd_calc_dir_to_target.tsv

or by default full atom rmsd with the internal lib:

	$ ./rna_calc_rmsd.py -t test_data/struc1.pdb -o test_output/rmsd_calc_dir_to_target.tsv test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb
	rmsd_calc_rmsd_to_target
	--------------------------------------------------------------------------------
	# of models: 4
	target:struc1.pdb      	rmsd_all
	struc1.pdb     	0.0
	struc2.pdb     	11.803
	struc3.pdb     	4.879
	struc4.pdb     	3.982
	tsv was created!  test_output/rmsd_calc_dir_to_target.tsv
	
.. select some residues and ignore some atoms 

    $ ./rna_calc_rmsd.py -t test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb --target_selection A:1-47+52-62 --model_selection A:1-47+52-62 /Users/magnus/work/src/rna-pdb-tools/rna_pdb_tools/utils/rna_calc_rmsd/test_data/pistol/clusters/pistol_thrs0.50A_clust01-000001_AA.pdb --model_ignore_selection A/57/O2\'  test_data/pistol/clusters/*_AA.pdb
    
    rmsd_calc_rmsd_to_target
    --------------------------------------------------------------------------------
     method:
    # of models: 4
    pistol_thrs0.50A_clust01-000001_AA.pdb 7.59648219355 1237
    pistol_thrs0.50A_clust01-000001_AA.pdb 7.59648219355 1237
    pistol_thrs0.50A_clust02-000001_AA.pdb 7.76647001671 1237
    pistol_thrs0.50A_clust03-000001_AA.pdb 18.1711015479 1237
    # of atoms used: 1237
    csv was created!  rmsds.csv

rmsd_calc_dir.py
-------------------------------------------------------------------------------

The program calculates all-atom rmsds (root-mean-square deviation) for all PDB structures in a given folder and save them to a file, as a matrix (which can be read by <https://github.com/mmagnus/rnastruc_clanstix>).

Usage:

    $ ./rna_calc_rmsd_all_vs_all.py -i test_data -o test_output/matrix.txt
    calc_rmsd_dir
    --------------------------------------------------------------------------------
     # of models: 4
    ... 1 test_data/struc1.pdb
    ... 2 test_data/struc2.pdb
    ... 3 test_data/struc3.pdb
    ... 4 test_data/struc4.pdb
    # test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb
    0.0 11.803 4.879 3.982
    11.803 0.0 12.153 12.043
    4.879 12.153 0.0 3.487
    3.982 12.043 3.487 0.0
    matrix was created!  test_output/matrix.txt

Output:

    $ head test_output/matrix.txt
    # test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb
    0.0 11.803 4.879 3.982
    11.803 0.0 12.153 12.043
    4.879 12.153 0.0 3.487
    3.982 12.043 3.487 0.0

The program is using Biopython:

Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878

The program is using (included in lib/) https://github.com/charnley/rmsd .
