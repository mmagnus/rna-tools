from Seq import RNASequence

seq = RNASequence("AAAUUAAGGGGAAGCGUUGAGCCGCUACCCAUAUGUGGUUCACUCGGAUAGCGGGGAGCUAAUAGUGAAACCGGCCCUUUAGGGG")
cst =             '...((((((((.(((......((((((.((....(((...)))..)).))))))...)))..............))))))))...'
seq.predict_ss("RNAfoldX", constraints=cst, verbose=True)
