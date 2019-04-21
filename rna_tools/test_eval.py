from Seq import RNASequence

seq = RNASequence("GGCAGGGGCGCUUCGGCCCCCUAUGCC")
seq.ss = "((((((((.((....)).)))).))))"
print(seq.eval())
