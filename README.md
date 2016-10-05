ClaRNA_play
-------------------------------------------------------------------------------

# Install 
Quick and dirty install (only ClaRNA will work)

To install ClaRNA_play, put add something like this in your `~/.bashrc` file

    export ClaRNAlib=/home/calmeida/Software/ClaRNA_play/ClaRNAlib/
    
ClaRNA requires:

    sudo pip install simplejson==2.6.1 networkx==1.8.1 scipy

# Start
To start:

    ./clarna_run.py -ipdb test_data/3e5fA_M15_c.pdb > test_output/3e5fA_M15_c.pdb.outCR
    ./clarna_run.py -ipdb test_data/3e5fA_M500_c.pdb > test_output/3e5fA_M500_c.pdb.outCR

to analyze:

    ./clarna_compare.py -iref test_output/3e5fA_M15_c.pdb.outCR -ichk test_output/3e5fA_M500_c.pdb.outCR
    3e5A_M15_c.pdb.outCR                       3e5fA_M500_c.pdb.outCR      0.662      NA         0.647      0.707      0.538      0.778      0.500      1.000

![](docs/clarna_run.png)

Tested with ClaRNA only.

    ./clarna_run.py -ipdb test_data/3e5fA_M500_c.pdb
    Classifier: Clarna
    A    5   A   49          bp C G                  WW_cis   0.6965
    A    6   A   48          bp C G                  WW_cis   0.6878
    A   10   A   34          bp A U                  WW_cis   0.6209
    A   10   A   35          bp A U                  WW_cis   0.6951
    A   13   A   28          bp A G                  WW_cis   0.9170
    A   14   A   27          bp U A                  WW_cis   0.7730
    A   15   A   26          bp G C                  WW_cis   0.8426
    A   16   A   25          bp G C                 WW_tran   0.6505
    A   17   A   24          bp C G                  WW_cis   0.7311
    A   18   A   23          bp G C                 WW_tran   0.7475
    A   31   A   36          bp G G                 HW_tran   0.6888
    A   41   A   46          bp C G                  WW_cis   0.8555

.. `clarna_compare.py -v`:
    
    ./clarna_compare.py -v -iref test_output/3e5fA_M15_c.pdb.outCR -ichk test_output/3e5fA_M500_c.pdb.outCR 
    test_output/3e5fA_M15_c.pdb.outCR
    test_output/3e5fA_M15_c.pdb.outCR
    A    2   A   52          bp U A                  WW_cis   0.9511
    A    3   A   51          bp U A                  WW_cis   0.8414
    A    4   A   50          bp C G                  WW_cis   0.7113
    A    5   A   49          bp C G                  WW_cis   0.9349
    A    6   A   48          bp C G                  WW_cis   0.8963
    A    9   A   35          bp A U                  WW_cis   0.7912
    A   10   A   34          bp A U                  WW_cis   0.9165
    A   11   A   33          bp G C                  WW_cis   0.8996
    A   13   A   27          bp A A                 WW_tran   0.6134
    A   13   A   28          bp A G                  WW_cis   0.8797
    A   15   A   26          bp G C                  WW_cis   0.8868
    A   16   A   25          bp G C                  WW_cis   0.7712
    A   17   A   24          bp C G                  WW_cis   0.9009
    A   18   A   23          bp G C                 WW_tran   0.7947
    A   19   A   22          bp G A                 SH_tran   0.7556
    A   31   A   36          bp G G                 HW_tran   0.8840
    A   38   A   48          bp A G                 WS_tran   0.6725
    A   40   A   47          bp C G                  WW_cis   0.9030
    A   41   A   46          bp C G                  WW_cis   0.8988

    test_output/3e5fA_M500_c.pdb.outCR
    A    5   A   49          bp C G                  WW_cis   0.6965
    A    6   A   48          bp C G                  WW_cis   0.6878
    A   10   A   34          bp A U                  WW_cis   0.6209
    A   10   A   35          bp A U                  WW_cis   0.6951
    A   13   A   28          bp A G                  WW_cis   0.9170
    A   14   A   27          bp U A                  WW_cis   0.7730
    A   15   A   26          bp G C                  WW_cis   0.8426
    A   16   A   25          bp G C                 WW_tran   0.6505
    A   17   A   24          bp C G                  WW_cis   0.7311
    A   18   A   23          bp G C                 WW_tran   0.7475
    A   31   A   36          bp G G                 HW_tran   0.6888
    A   41   A   46          bp C G                  WW_cis   0.8555

    compare:
    3e5fA_M15_c.pdb.outCR     3e5fA_M500_c.pdb.outCR

    reference vs prediction:                                                                                                                                                                                                                                                       
      WC    ref    13   pred-match     7                                                                                                                                                                                                                                          
      nWC    ref     6   pred-match     3                                                                                                                                                                                                                                          
    stack    ref     0   pred-match     0                                                                                                                                                                                                                                          
																																
    prediction vs reference:   WC   pred     8    ref-match     6                                                                                                                                                                                                                  
      nWC   pred     4    ref-match     4                                                                                                                                                                                                                                          
    stack   pred     0    ref-match     0                                                                                                                                                                                                                                          
																																
    for all motifs:nref   19,  nchkTP   10,  nFN    9                                                                                                                                                                                                                              
    nchk   12,  nrefTP   10,  nFP    2                                                                                                                                                                                                                                             
																																
																																
    summary:                                                                                                                                                                                                                                                                       
    inf_all      0.662                                                                                                                                                                                                                                                             
    inf_stack -999.999                                                                                                                                                                                                                                                             
    inf_WC       0.647                                                                                                                                                                                                                                                             
    inf_nWC      0.707                                                                                                                                                                                                                                                             
    SNS_WC       0.538                                                                                                                                                                                                                                                             
    PPV_WC       0.778                                                                                                                                                                                                                                                             
    SNS_nWC      0.500                                                                                                                                                                                                                                                             
    PPV_nWC      1.000