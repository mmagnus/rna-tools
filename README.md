ClaRNA_play
-------------------------------------------------------------------------------

    clarna_run.py -ipdb 1xjr.pdb > 1xjr.pdb.CL

    clarna_anal.py -iref 1xjr.pdb.outCR -ichk 1xjr.pdb.outCR
    1xjr.pdb.outCR                               1xjr.pdb.outCR      1.000      NA         1.000      1.000      1.000      1.000      1.000      1.000


![](docs/clarna_run.png)

Tested with ClaRNA only.

    ➜  magnus clarna_run.py -ipdb 1xjr.pdb
    Classifier: Clarna
    A    4   A   46          bp G U                  WW_cis   0.7957
    A    5   A   45          bp U A                  WW_cis   0.9186
    A    6   A   44          bp U A                  WW_cis   0.9299
    A    7   A   43          bp C G                  WW_cis   0.9177
    A    8   A   42          bp A U                  WW_cis   0.9332
    A    9   A   41          bp C G                  WW_cis   0.9209
    A   10   A   40          bp C A                  WW_cis   0.6129
    A   11   A   33          bp G A                  SS_cis   0.7598
    A   13   A   39          bp G C                  WW_cis   0.8812
    A   14   A   37          bp G U                  WW_cis   0.8378
    A   15   A   36          bp C G                  WW_cis   0.9079
    A   16   A   35          bp C G                  WW_cis   0.9229
    A   17   A   34          bp A G                  WW_cis   0.8822
    A   18   A   32          bp C G                  WW_cis   0.9262
    A   19   A   31          bp G C                  WW_cis   0.8962
    A   19   A   20          bp G C                  SS_cis   0.7453
    A   20   A   28          bp C G                  WW_cis   0.7902
    A   21   A   27          bp G C                  WW_cis   0.8676
    A   22   A   26          bp G A                 SH_tran   0.7988
    A   28   A   31          bp G C                  SS_cis   0.6844
    A   38   A   39          bp A C                  SH_cis   0.8775

