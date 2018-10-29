RNA Secondary Structure prediction using MC-Fold including mask (constraints)
-------------------------------------------------------------------------------

181029

The code to run this analysis is:

    cat test_seq_mcfold.py
    #!/usr/bin/python
    #-*- coding: utf-8 -*-
    from rna_pdb_tools.Seq import RNASequence

    seq = RNASequence("CCuuuuGG")
    print(seq.predict_ss("mcfold", verbose=True))

    seq = RNASequence("CCuuuuGG")
    print(seq.predict_ss("mcfold", constraints='(......)', verbose=True))

to run it:

    ./test_seq_mcfold.py
    curl -Y 0 -y 300 -F "pass=lucy" -F sequence="CCuuuuGG" http://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
    (-5.43, '(((..)))')
    curl -Y 0 -y 300 -F "pass=lucy" -F mask="(......)" -F sequence="CCuuuuGG" http://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
    (0.0, '')

and from the console if you want:

    [mm] tests$ git:(master) ✗ ./test_rna_ss_pred.sh
    + rna_ss_pred.py --method mcfold --seq gggaaacc
    gggaaacc
    (((..))) -7.61
    + rna_ss_pred.py --method mcfold --seq gggaaacc --cst '((....))'
    gggaaacc
    ((....)) <= cst
    x 0.0
    + rna_ss_pred.py --method mcfold --seq gggaaacc --cst '(((..)))'
    gggaaacc
    (((..))) <= cst
    (((..))) -7.61
    + rna_ss_pred.py --method mcfold --file rna1.fa --cstinfile
    gggaaacc
    x 0.0
    + rna_ss_pred.py --method mcfold --file rna1.fa
    gggaaacc
    (((..))) -7.61

the distinction between unpaired and not constrained:

    [mm] Desktop$ rna_ss_pred.py --method mcfold --file seqs/giic+g.fa --cstinfile
    gGUAUguaaGUACc
    (((((....))))) <= cst
    (((((....))))) -12.49
    [mm] Desktop$ rna_ss_pred.py --method mcfold --file seqs/giic+g.fa --cstinfile
    gGUAUguaaGUACc
    (((((****))))) <= cst
    ((((((..)))))) -13.62
    
and an interesting issue of GC vs CG at the end:

    [mm] Desktop$ rna_ss_pred.py --method mcfold --file seqs/upper4x.fa
    gCCCUguaaAGGGc
    ((((((..)))))) -16.07
    [mm] Desktop$ rna_ss_pred.py --method mcfold --file seqs/upper4x.fa
    cCCCUguaaAGGGg
    ((((((..)))))) -15.74

