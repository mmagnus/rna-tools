rna_ss_pred.py
-------------------------------------------------------------------------------

181026

`rna_ss_pred.py` is a cmd tool that uses wrappers implemented in rna-pdb-tools.

    rna_ss_pred.py --cyclefold --seq GUAUguaaGUAC
    GUAUguaaGUAC
    (((((..))))) -9.9541

    rna_ss_pred.py --mcfold    --seq GUAUguaaGUAC
    (((((..)))))  -9.97

    rna_ss_pred.py --cyclefold --seq UGGCguaaGACA
    UGGCguaaGACA
    (((((..))))) -10.842

    rna_ss_pred.py --mcfold    --seq UGGCguaaGACA
    (((((..))))) -10.82
