from Seq import RNASequence

# main
if __name__ == '__main__':
    seq = 'acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg'
    ss =  '((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))'
    seq = RNASequence(seq)
    seq.ss = ss
    fe = seq.eval()
    print(fe)
    fa = seq.get_foldability()
    print(fa)

    seq = RNASequence("GGCAGGGGCGCUUCGGCCCCCUAUGCC")
    seq.ss =          "((((((((.((....)).)))).))))"
    print(seq.eval())

    fa = seq.get_foldability()
    print(fa)
