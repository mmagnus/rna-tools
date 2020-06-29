rna_calc_inf
-------------------------------------------------------------------------------

````
rna_calc_inf.py -t ../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb \
                ../rna_calc_rmsd/test_data/pistol/clusters/*.pdb \
                -o test_data/pistol_inf.csv --print-results
100% (4 of 4) |#########################################################################################################################################| Elapsed Time: 0:00:00 ETA:  00:00:00
csv was created!  test_data/pistol_inf.csv
                                           target                                              fn  inf_all  inf_stack  inf_WC  inf_nWC  sns_WC  ppv_WC  sns_nWC  ppv_nWC
0  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb          pistol_thrs0.50A_clust02-000001_AA.pdb     0.41        0.0    0.86     0.38    0.89    0.84     0.33     0.43
1  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb          pistol_thrs0.50A_clust01-000001_AA.pdb     0.42        0.0    0.87     0.35    0.94    0.81     0.33     0.38
2  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb     0.57        0.0    1.00     1.00    1.00    1.00     1.00     1.00
3  5k7c_clean_onechain_renumber_as_puzzle_srr.pdb          pistol_thrs0.50A_clust03-000001_AA.pdb     0.37        0.0    0.75     0.30    0.83    0.68     0.22     0.40
```

```
rna_calc_inf.py -t 17_0_solution_5K7C_rpr.pdb --target-selection A:1-47+52-62 --model-selection A:1-47+52-62 *.pdb --print-results --sort-results

100% (109 of 109) |####################################################################################################################################################################| Elapsed Time: 0:17:07 ETA:  00:00:00
csv was created!  inf.csv
                             target                                             fn  inf_all  inf_stack  inf_WC  inf_nWC  sns_WC  ppv_WC  sns_nWC  ppv_nWC
2    17_0_solution_5K7C_rpr_sel.pdb  17_0_solution_5K7C_MissAtomResi53_rpr_sel.pdb     1.00       1.00    1.00     1.00    1.00    1.00     1.00     1.00
5    17_0_solution_5K7C_rpr_sel.pdb                 17_0_solution_5K7C_rpr_sel.pdb     1.00       1.00    1.00     1.00    1.00    1.00     1.00     1.00
27   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_1_rpr_sel.pdb     0.81       0.83    0.92     0.41    0.94    0.90     0.33     0.50
23   17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_2_rpr_sel.pdb     0.77       0.79    0.84     0.38    0.89    0.80     0.22     0.67
37   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_4_rpr_sel.pdb     0.76       0.78    0.86     0.33    0.89    0.84     0.22     0.50
17   17_0_solution_5K7C_rpr_sel.pdb                  17_DasExtraInfo_1_rpr_sel.pdb     0.76       0.78    0.84     0.41    0.89    0.80     0.33     0.50
49   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_8_rpr_sel.pdb     0.76       0.77    0.90     0.30    0.94    0.85     0.22     0.40
20   17_0_solution_5K7C_rpr_sel.pdb                  17_DasExtraInfo_2_rpr_sel.pdb     0.76       0.77    0.86     0.38    0.89    0.84     0.22     0.67
46   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_7_rpr_sel.pdb     0.75       0.78    0.83     0.33    0.83    0.83     0.22     0.50
52   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_9_rpr_sel.pdb     0.75       0.78    0.86     0.15    0.89    0.84     0.11     0.20
43   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_6_rpr_sel.pdb     0.75       0.77    0.84     0.38    0.89    0.80     0.22     0.67
3    17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_6_rpr_sel.pdb     0.75       0.76    0.84     0.47    0.89    0.80     0.22     1.00
9    17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_8_rpr_sel.pdb     0.75       0.78    0.89     0.17    0.89    0.89     0.11     0.25
89   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_1_rpr_sel.pdb     0.75       0.78    0.93     0.13    1.00    0.86     0.11     0.14
6    17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_7_rpr_sel.pdb     0.74       0.75    0.89     0.19    0.89    0.89     0.11     0.33
40   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_5_rpr_sel.pdb     0.74       0.77    0.84     0.17    0.89    0.80     0.11     0.25
22   17_0_solution_5K7C_rpr_sel.pdb                  17_DasExtraInfo_3_rpr_sel.pdb     0.74       0.77    0.84     0.24    0.89    0.80     0.11     0.50
103  17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_7_rpr_sel.pdb     0.73       0.76    0.84     0.33    0.89    0.80     0.33     0.33
30   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_2_rpr_sel.pdb     0.73       0.76    0.86     0.15    0.89    0.84     0.11     0.20
105  17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_8_rpr_sel.pdb     0.72       0.74    0.82     0.38    0.89    0.76     0.33     0.43
34   17_0_solution_5K7C_rpr_sel.pdb                           17_Das_3_rpr_sel.pdb     0.72       0.74    0.83     0.38    0.83    0.83     0.33     0.43
4    17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_6_rpr_sel.pdb     0.72       0.74    0.73     0.58    0.67    0.80     0.33     1.00
93   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_2_rpr_sel.pdb     0.71       0.73    0.77     0.38    0.89    0.67     0.22     0.67
14   17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_9_rpr_sel.pdb     0.71       0.72    0.86     0.19    0.83    0.88     0.11     0.33
28   17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_4_rpr_sel.pdb     0.71       0.72    0.86     0.17    0.89    0.84     0.11     0.25
7    17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_7_rpr_sel.pdb     0.70       0.71    0.83     0.24    0.83    0.83     0.11     0.50
29   17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_4_rpr_sel.pdb     0.70       0.72    0.76     0.50    0.67    0.86     0.33     0.75
19   17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_1_rpr_sel.pdb     0.70       0.72    0.80     0.19    0.78    0.82     0.11     0.33
26   17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_3_rpr_sel.pdb     0.70       0.73    0.81     0.17    0.83    0.79     0.11     0.25
101  17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_6_rpr_sel.pdb     0.70       0.74    0.78     0.24    0.78    0.78     0.22     0.25
97   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_4_rpr_sel.pdb     0.70       0.74    0.74     0.33    0.78    0.70     0.22     0.50
68   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_2_rpr_sel.pdb     0.70       0.77    0.80     0.20    0.78    0.82     0.22     0.18
107  17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_9_rpr_sel.pdb     0.70       0.71    0.82     0.19    0.89    0.76     0.11     0.33
94   17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_3_rpr_sel.pdb     0.69       0.75    0.75     0.00    0.61    0.92     0.00     0.00
81   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_6_rpr_sel.pdb     0.69       0.75    0.77     0.21    0.72    0.81     0.22     0.20
95   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_3_rpr_sel.pdb     0.69       0.70    0.78     0.41    0.78    0.78     0.33     0.50
71   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_3_rpr_sel.pdb     0.69       0.75    0.74     0.22    0.72    0.76     0.22     0.22
99   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS2_5_rpr_sel.pdb     0.69       0.72    0.77     0.27    0.72    0.81     0.22     0.33
0    17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_5_rpr_sel.pdb     0.69       0.72    0.81     0.15    0.83    0.79     0.11     0.20
21   17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_2_rpr_sel.pdb     0.69       0.69    0.84     0.17    0.89    0.80     0.11     0.25
32   17_0_solution_5K7C_rpr_sel.pdb                     17_Dohkolyan_3_rpr_sel.pdb     0.69       0.71    0.83     0.00    0.83    0.83     0.00     0.00
13   17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_9_rpr_sel.pdb     0.69       0.72    0.82     0.15    0.89    0.76     0.11     0.20
15   17_0_solution_5K7C_rpr_sel.pdb                     17_Bujnicki_10_rpr_sel.pdb     0.69       0.69    0.84     0.17    0.89    0.80     0.11     0.25
18   17_0_solution_5K7C_rpr_sel.pdb                      17_Bujnicki_1_rpr_sel.pdb     0.69       0.72    0.82     0.15    0.89    0.76     0.11     0.20
25   17_0_solution_5K7C_rpr_sel.pdb                          17_Das_10_rpr_sel.pdb     0.69       0.68    0.84     0.33    0.89    0.80     0.22     0.50
88   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_8_rpr_sel.pdb     0.68       0.72    0.74     0.25    0.72    0.76     0.22     0.29
58   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_1_rpr_sel.pdb     0.68       0.77    0.59     0.19    0.56    0.62     0.11     0.33
1    17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_5_rpr_sel.pdb     0.68       0.69    0.72     0.50    0.61    0.85     0.33     0.75
39   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_5_rpr_sel.pdb     0.68       0.74    0.67     0.19    0.67    0.67     0.11     0.33
45   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_7_rpr_sel.pdb     0.68       0.71    0.77     0.19    0.72    0.81     0.11     0.33
36   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_4_rpr_sel.pdb     0.68       0.75    0.63     0.33    0.61    0.65     0.22     0.50
65   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_1_rpr_sel.pdb     0.67       0.71    0.74     0.25    0.72    0.76     0.22     0.29
10   17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_8_rpr_sel.pdb     0.67       0.66    0.86     0.00    0.83    0.88     0.00     0.00
48   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_8_rpr_sel.pdb     0.67       0.66    0.80     0.33    0.78    0.82     0.22     0.50
92   17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_2_rpr_sel.pdb     0.67       0.68    0.76     0.33    0.67    0.86     0.22     0.50
70   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_2_rpr_sel.pdb     0.66       0.64    0.86     0.24    0.83    0.88     0.11     0.50
77   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_5_rpr_sel.pdb     0.66       0.72    0.71     0.20    0.67    0.75     0.22     0.18
24   17_0_solution_5K7C_rpr_sel.pdb                          17_Chen_3_rpr_sel.pdb     0.65       0.68    0.73     0.00    0.67    0.80     0.00     0.00
55   17_0_solution_5K7C_rpr_sel.pdb                         17_Ding_10_rpr_sel.pdb     0.64       0.65    0.77     0.30    0.72    0.81     0.22     0.40
51   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_9_rpr_sel.pdb     0.64       0.65    0.77     0.17    0.72    0.81     0.11     0.25
74   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_4_rpr_sel.pdb     0.64       0.71    0.71     0.19    0.67    0.75     0.22     0.17
12   17_0_solution_5K7C_rpr_sel.pdb                       17_Adamiak_3_rpr_sel.pdb     0.64       0.61    0.81     0.38    0.83    0.79     0.22     0.67
33   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_3_rpr_sel.pdb     0.64       0.69    0.63     0.19    0.61    0.65     0.11     0.33
98   17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_5_rpr_sel.pdb     0.63       0.66    0.71     0.00    0.56    0.91     0.00     0.00
82   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_6_rpr_sel.pdb     0.62       0.61    0.76     0.24    0.78    0.74     0.11     0.50
42   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_6_rpr_sel.pdb     0.62       0.64    0.73     0.17    0.67    0.80     0.11     0.25
85   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_7_rpr_sel.pdb     0.62       0.63    0.76     0.19    0.78    0.74     0.11     0.33
38   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_1_rpr_sel.pdb     0.62       0.66    0.65     0.25    0.61    0.69     0.22     0.29
76   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_4_rpr_sel.pdb     0.61       0.60    0.76     0.24    0.78    0.74     0.11     0.50
79   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_5_rpr_sel.pdb     0.61       0.58    0.82     0.24    0.78    0.88     0.11     0.50
8    17_0_solution_5K7C_rpr_sel.pdb                       17_Adamiak_1_rpr_sel.pdb     0.61       0.57    0.79     0.38    0.83    0.75     0.22     0.67
11   17_0_solution_5K7C_rpr_sel.pdb                       17_Adamiak_2_rpr_sel.pdb     0.60       0.59    0.76     0.24    0.78    0.74     0.11     0.50
57   17_0_solution_5K7C_rpr_sel.pdb                     17_Dohkolyan_2_rpr_sel.pdb     0.60       0.62    0.73     0.17    0.67    0.80     0.11     0.25
87   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_8_rpr_sel.pdb     0.59       0.57    0.77     0.24    0.72    0.81     0.11     0.50
66   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_1_rpr_sel.pdb     0.59       0.66    0.50     0.27    0.44    0.57     0.22     0.33
47   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_4_rpr_sel.pdb     0.59       0.66    0.61     0.00    0.61    0.61     0.00     0.00
31   17_0_solution_5K7C_rpr_sel.pdb                          17_Ding_2_rpr_sel.pdb     0.59       0.62    0.63     0.24    0.61    0.65     0.11     0.50
90   17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_1_rpr_sel.pdb     0.59       0.60    0.75     0.00    0.61    0.92     0.00     0.00
80   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_6_rpr_sel.pdb     0.58       0.64    0.50     0.33    0.44    0.57     0.22     0.50
50   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_5_rpr_sel.pdb     0.58       0.64    0.57     0.25    0.56    0.59     0.22     0.29
44   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_3_rpr_sel.pdb     0.58       0.64    0.61     0.13    0.61    0.61     0.11     0.14
106  17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_9_rpr_sel.pdb     0.58       0.61    0.63     0.24    0.44    0.89     0.11     0.50
63   17_0_solution_5K7C_rpr_sel.pdb               17_RNAComposerAS1_10_rpr_sel.pdb     0.58       0.65    0.50     0.27    0.44    0.57     0.22     0.33
64   17_0_solution_5K7C_rpr_sel.pdb               17_RNAComposerAS2_10_rpr_sel.pdb     0.58       0.58    0.68     0.24    0.56    0.83     0.11     0.50
67   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_1_rpr_sel.pdb     0.58       0.54    0.79     0.24    0.83    0.75     0.11     0.50
104  17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_8_rpr_sel.pdb     0.57       0.61    0.58     0.00    0.33    1.00     0.00     0.00
108  17_0_solution_5K7C_rpr_sel.pdb                         17_Xiao_10_rpr_sel.pdb     0.57       0.61    0.58     0.00    0.39    0.88     0.00     0.00
61   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_9_rpr_sel.pdb     0.56       0.58    0.55     0.45    0.50    0.60     0.33     0.60
75   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_4_rpr_sel.pdb     0.56       0.61    0.50     0.33    0.44    0.57     0.22     0.50
41   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_2_rpr_sel.pdb     0.56       0.61    0.60     0.15    0.61    0.58     0.11     0.20
72   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_3_rpr_sel.pdb     0.56       0.60    0.59     0.24    0.50    0.69     0.22     0.25
54   17_0_solution_5K7C_rpr_sel.pdb                     17_Dohkolyan_1_rpr_sel.pdb     0.56       0.64    0.50     0.14    0.44    0.57     0.11     0.17
86   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_8_rpr_sel.pdb     0.55       0.60    0.50     0.30    0.44    0.57     0.22     0.40
62   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_9_rpr_sel.pdb     0.55       0.54    0.72     0.19    0.61    0.85     0.11     0.33
73   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS2_3_rpr_sel.pdb     0.55       0.59    0.55     0.30    0.39    0.78     0.22     0.40
96   17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_4_rpr_sel.pdb     0.55       0.60    0.58     0.17    0.39    0.88     0.11     0.25
102  17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_7_rpr_sel.pdb     0.53       0.55    0.62     0.00    0.39    1.00     0.00     0.00
91   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_9_rpr_sel.pdb     0.53       0.64    0.30     0.27    0.28    0.31     0.22     0.33
78   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_5_rpr_sel.pdb     0.53       0.56    0.59     0.25    0.50    0.69     0.22     0.29
100  17_0_solution_5K7C_rpr_sel.pdb                          17_Xiao_6_rpr_sel.pdb     0.53       0.55    0.62     0.00    0.39    1.00     0.00     0.00
59   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_8_rpr_sel.pdb     0.50       0.56    0.47     0.17    0.33    0.67     0.11     0.25
53   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_6_rpr_sel.pdb     0.49       0.58    0.43     0.00    0.33    0.55     0.00     0.00
84   17_0_solution_5K7C_rpr_sel.pdb                     17_SimRNAAS1_7_rpr_sel.pdb     0.48       0.62    0.17     0.22    0.17    0.18     0.22     0.22
69   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_2_rpr_sel.pdb     0.48       0.60    0.20     0.21    0.17    0.25     0.22     0.20
56   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_7_rpr_sel.pdb     0.48       0.55    0.44     0.00    0.28    0.71     0.00     0.00
83   17_0_solution_5K7C_rpr_sel.pdb                17_RNAComposerAS1_7_rpr_sel.pdb     0.48       0.61    0.20     0.20    0.17    0.23     0.22     0.18
60   17_0_solution_5K7C_rpr_sel.pdb                         17_Major_9_rpr_sel.pdb     0.45       0.50    0.42     0.17    0.22    0.80     0.11     0.25
16   17_0_solution_5K7C_rpr_sel.pdb                         17_Chen_10_rpr_sel.pdb     0.45       0.56    0.20     0.21    0.17    0.23     0.22     0.20
35   17_0_solution_5K7C_rpr_sel.pdb                        17_Major_10_rpr_sel.pdb     0.41       0.48    0.36     0.14    0.22    0.57     0.11     0.17
```
