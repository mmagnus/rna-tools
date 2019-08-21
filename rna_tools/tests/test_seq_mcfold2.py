from rna_tools.Seq import RNASequence

seq = "AGUCGUUGUGGCGACUAUAACCAAGCUCUUUAAGCCACAAGCGUUGCUGAUGAGGUUUCAUAACAUCAGCAGGUAGAG"
ss =  "((((..[[[[[.)))).........((((.....]]]]]....(((((((((...........)))))))))..))))"
s = RNASequence(seq, ss)
print(s.predict_ss(method='mcfold', constraints=ss, verbose=True))
ss =  "((((..[[[[[.)))).........((((.....]]]]]....(((((((((***********)))))))))..))))"
print(s.predict_ss(method='mcfold', constraints=ss, verbose=True))
