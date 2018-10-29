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
