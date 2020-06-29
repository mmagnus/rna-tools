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



```
rna_calc_rmsd_multi_targets.py --models *.pdb \
                                    --targets solutions/*.pdb \
                                    --target-selection A:1-27+29-41 \
                                    --model-selection A:1-27+29-41
solutions/21_solution_0_ChainA.pdb
--------------------------------------------------------------------------------
# method: all-atom-built-in
# of models: 60
# of atoms used: 861
csv was created!  21_solution_0_ChainA.csv
solutions/21_solution_0_ChainB.pdb
--------------------------------------------------------------------------------
# method: all-atom-built-in
# of models: 60
# of atoms used: 861
csv was created!  21_solution_0_ChainB.csv
solutions/21_solution_1_ChainA.pdb
--------------------------------------------------------------------------------
# method: all-atom-built-in
# of models: 60
# of atoms used: 861
csv was created!  21_solution_1_ChainA.csv
solutions/21_solution_1_ChainB.pdb
--------------------------------------------------------------------------------
# method: all-atom-built-in
# of models: 60
# of atoms used: 861
csv was created!  21_solution_1_ChainB.csv
solutions/21_solution_2.pdb
--------------------------------------------------------------------------------
# method: all-atom-built-in
# of models: 60
# of atoms used: 861
csv was created!  21_solution_2.csv
--------------------------------------------------------------------------------
                          21_solution_0_ChainA.pdb  21_solution_0_ChainB.pdb  21_solution_1_ChainA.pdb  21_solution_1_ChainB.pdb  21_solution_2.pdb   mean  median  variance    min    max    sd
fn
21_3dRNA_1_rpr.pdb                           12.17                     12.11                     12.17                     12.11              12.11  12.13   12.11      0.00  12.11  12.17  0.03
21_3dRNA_2_rpr.pdb                           12.48                     12.43                     12.48                     12.43              12.40  12.44   12.43      0.00  12.40  12.48  0.04
21_3dRNA_3_rpr.pdb                           12.04                     12.05                     12.04                     12.05              12.13  12.06   12.05      0.00  12.04  12.13  0.04
21_3dRNA_4_rpr.pdb                           18.01                     17.96                     18.01                     17.96              18.06  18.00   18.01      0.00  17.96  18.06  0.04
21_3dRNA_5_rpr.pdb                           17.61                     17.49                     17.61                     17.49              17.57  17.55   17.57      0.00  17.49  17.61  0.06
21_Adamiak_1_rpr.pdb                          4.64                      4.61                      4.64                      4.61               4.64   4.63    4.64      0.00   4.61   4.64  0.02
21_Adamiak_2_rpr.pdb                          4.65                      4.62                      4.65                      4.62               4.64   4.64    4.64      0.00   4.62   4.65  0.02
21_Adamiak_3_rpr.pdb                          4.61                      4.58                      4.61                      4.58               4.58   4.59    4.58      0.00   4.58   4.61  0.02
21_Adamiak_4_rpr.pdb                          4.58                      4.51                      4.58                      4.51               4.54   4.54    4.54      0.00   4.51   4.58  0.04
21_Adamiak_5_rpr.pdb                          5.50                      5.56                      5.50                      5.56               5.56   5.54    5.56      0.00   5.50   5.56  0.03
21_Bujnicki_1_rpr.pdb                         6.77                      6.77                      6.77                      6.77               6.78   6.77    6.77      0.00   6.77   6.78  0.00
21_Bujnicki_2_rpr.pdb                         8.69                      8.84                      8.69                      8.84               8.79   8.77    8.79      0.01   8.69   8.84  0.08
21_Bujnicki_3_rpr.pdb                        11.98                     11.97                     11.98                     11.97              12.20  12.02   11.98      0.01  11.97  12.20  0.10
21_Bujnicki_4_rpr.pdb                         6.44                      6.39                      6.44                      6.39               6.41   6.41    6.41      0.00   6.39   6.44  0.03
21_Bujnicki_5_rpr.pdb                        21.56                     21.54                     21.56                     21.54              21.69  21.58   21.56      0.00  21.54  21.69  0.06
21_ChenHighLig_1_rpr.pdb                      4.01                      3.97                      4.01                      3.97               4.07   4.01    4.01      0.00   3.97   4.07  0.04
21_ChenHighLig_2_rpr.pdb                      3.79                      3.75                      3.79                      3.75               3.71   3.76    3.75      0.00   3.71   3.79  0.03
21_ChenHighLig_3_rpr.pdb                      5.08                      4.95                      5.08                      4.95               5.01   5.01    5.01      0.00   4.95   5.08  0.07
21_ChenHighLig_4_rpr.pdb                      4.59                      4.59                      4.59                      4.59               4.65   4.60    4.59      0.00   4.59   4.65  0.03
21_ChenHighLig_5_rpr.pdb                      7.00                      6.97                      7.00                      6.97               6.96   6.98    6.97      0.00   6.96   7.00  0.02
21_ChenLowLig_1_rpr.pdb                       4.01                      3.97                      4.01                      3.97               4.07   4.01    4.01      0.00   3.97   4.07  0.04
21_ChenLowLig_2_rpr.pdb                       3.79                      3.75                      3.79                      3.75               3.71   3.76    3.75      0.00   3.71   3.79  0.03
21_ChenLowLig_3_rpr.pdb                       5.08                      4.95                      5.08                      4.95               5.01   5.01    5.01      0.00   4.95   5.08  0.07
21_ChenLowLig_4_rpr.pdb                       4.59                      4.59                      4.59                      4.59               4.65   4.60    4.59      0.00   4.59   4.65  0.03
21_ChenLowLig_5_rpr.pdb                       7.00                      6.97                      7.00                      6.97               6.96   6.98    6.97      0.00   6.96   7.00  0.02
21_DasLORES_1_rpr.pdb                         4.24                      4.10                      4.24                      4.10               4.10   4.16    4.10      0.01   4.10   4.24  0.08
21_DasLORES_2_rpr.pdb                         4.51                      4.45                      4.51                      4.45               4.49   4.48    4.49      0.00   4.45   4.51  0.03
21_DasLORES_3_rpr.pdb                        14.46                     14.26                     14.46                     14.26              14.42  14.37   14.42      0.01  14.26  14.46  0.10
21_DasLORES_4_rpr.pdb                         6.39                      6.43                      6.39                      6.43               6.50   6.43    6.43      0.00   6.39   6.50  0.04
21_DasLORES_5_rpr.pdb                         8.92                      8.98                      8.92                      8.98               9.16   8.99    8.98      0.01   8.92   9.16  0.10
21_Das_1_rpr.pdb                              5.71                      5.60                      5.71                      5.60               5.61   5.65    5.61      0.00   5.60   5.71  0.06
21_Das_2_rpr.pdb                              5.53                      5.42                      5.53                      5.42               5.44   5.47    5.44      0.00   5.42   5.53  0.06
21_Das_3_rpr.pdb                              5.07                      5.06                      5.07                      5.06               5.01   5.05    5.06      0.00   5.01   5.07  0.03
21_Das_4_rpr.pdb                              5.73                      5.57                      5.73                      5.57               5.61   5.64    5.61      0.01   5.57   5.73  0.08
21_Das_5_rpr.pdb                              6.04                      6.01                      6.04                      6.01               6.13   6.05    6.04      0.00   6.01   6.13  0.05
21_RNAComposer_1_rpr.pdb                     14.79                     14.70                     14.79                     14.70              14.69  14.73   14.70      0.00  14.69  14.79  0.05
21_RNAComposer_2_rpr.pdb                     18.32                     18.14                     18.32                     18.14              18.24  18.23   18.24      0.01  18.14  18.32  0.09
21_RNAComposer_3_rpr.pdb                     12.30                     12.12                     12.30                     12.12              12.26  12.22   12.26      0.01  12.12  12.30  0.09
21_RNAComposer_4_rpr.pdb                     13.24                     12.99                     13.24                     12.99              13.13  13.12   13.13      0.02  12.99  13.24  0.13
21_RNAComposer_5_rpr.pdb                     17.03                     16.94                     17.03                     16.94              17.02  16.99   17.02      0.00  16.94  17.03  0.05
21_RW3D_1_rpr.pdb                            18.10                     17.85                     18.10                     17.85              18.01  17.98   18.01      0.02  17.85  18.10  0.13
21_RW3D_2_rpr.pdb                            16.29                     16.06                     16.29                     16.06              16.19  16.18   16.19      0.01  16.06  16.29  0.12
21_RW3D_3_rpr.pdb                            17.91                     17.73                     17.91                     17.73              17.87  17.83   17.87      0.01  17.73  17.91  0.09
21_RW3D_4_rpr.pdb                            18.38                     18.14                     18.38                     18.14              18.29  18.27   18.29      0.01  18.14  18.38  0.12
21_RW3D_5_rpr.pdb                            18.29                     18.05                     18.29                     18.05              18.21  18.18   18.21      0.01  18.05  18.29  0.12
21_RW3D_6_rpr.pdb                            16.98                     16.86                     16.98                     16.86              16.96  16.93   16.96      0.00  16.86  16.98  0.06
21_RW3D_7_rpr.pdb                            18.27                     18.01                     18.27                     18.01              18.17  18.15   18.17      0.02  18.01  18.27  0.13
21_RW3D_8_rpr.pdb                            17.83                     17.60                     17.83                     17.60              17.75  17.72   17.75      0.01  17.60  17.83  0.12
21_RW3D_9_rpr.pdb                            18.15                     17.93                     18.15                     17.93              18.08  18.05   18.08      0.01  17.93  18.15  0.11
21_RW3D_10_rpr.pdb                           16.20                     16.08                     16.20                     16.08              16.18  16.15   16.18      0.00  16.08  16.20  0.06
21_Sanbonmatsu_1_rpr.pdb                     14.84                     14.61                     14.84                     14.61              14.73  14.73   14.73      0.01  14.61  14.84  0.12
21_Sanbonmatsu_2_rpr.pdb                     14.84                     14.61                     14.84                     14.61              14.73  14.73   14.73      0.01  14.61  14.84  0.12
21_Sanbonmatsu_3_rpr.pdb                     13.95                     14.11                     13.95                     14.11              14.13  14.05   14.11      0.01  13.95  14.13  0.09
21_Sanbonmatsu_4_rpr.pdb                     13.95                     14.11                     13.95                     14.11              14.13  14.05   14.11      0.01  13.95  14.13  0.09
21_simRNA_1_rpr.pdb                          11.88                     11.88                     11.88                     11.88              12.10  11.92   11.88      0.01  11.88  12.10  0.10
21_simRNA_2_rpr.pdb                          11.58                     11.55                     11.58                     11.55              11.71  11.59   11.58      0.00  11.55  11.71  0.07
21_simRNA_3_rpr.pdb                          13.15                     13.15                     13.15                     13.15              13.36  13.19   13.15      0.01  13.15  13.36  0.09
21_simRNA_4_rpr.pdb                          10.65                     10.65                     10.65                     10.65              10.75  10.67   10.65      0.00  10.65  10.75  0.04
21_simRNA_5_rpr.pdb                          13.15                     13.15                     13.15                     13.15              13.36  13.19   13.15      0.01  13.15  13.36  0.09
21_solution_2_rpr.pdb                         0.87                      0.90                      0.87                      0.90               0.00   0.71    0.87      0.16   0.00   0.90  0.40
--------------------------------------------------------------------------------
Save rna_calc_rmsd_multi_targets_output.csv
```
