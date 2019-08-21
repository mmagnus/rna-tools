#!/usr/bin/python
#-*- coding: utf-8 -*-
from rna_tools.Seq import RNASequence

seq = RNASequence("CCuuuuGG")
print(seq.predict_ss("mcfold", verbose=True))

seq = RNASequence("CCuuuuGG")
print(seq.predict_ss("mcfold", constraints='(......)', verbose=True))

import rna_tools.tools.rna_alignment.rna_alignment as ra

alignment = ra.RNAalignment('../input/RF02679_pistol_stockholm.txt')
for s in alignment:
    s.remove_gaps()
    print(s)
    print(s.seq)
    print(s.ss)
    # hack
    s2 = RNASequence(s.seq, s.ss)
    ss = s2.predict_ss(method='mcfold', constraints=s.ss.replace('.', '*'), verbose=True)
    print(ss)
    break

seq = 'acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg'
ss =  '((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))'
s = RNASequence(seq, ss)
print(s)

cst1 = "((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))"
cst2 = "((((........)))).......((((..............((((((((((..............))))))))))..))))"
cst3 = "((((xxxxxxxx))))xxxxxxx((((xxxxxxxxxxxxxx((((((((((xxxxxxxxxxxxxx))))))))))xx))))"
cst4 = "((((********))))*******((((**************((((((((((**************))))))))))**))))"
cst5 = "((((..[[[[[..)))).......((((....]]]]]....(((((((((((((((....)))))))))))))))..))))".replace('.', '*')
## ss = s.predict_ss(method='mcfold', constraints=cst1, verbose=True)
## print(ss)
## ss = s.predict_ss(method='mcfold', constraints=cst2, verbose=True)
## print(ss)
## ss = s.predict_ss(method='mcfold', constraints=cst3, verbose=True)
## print(ss)
for cst in [cst5]:
    ss = s.predict_ss(method='mcfold', constraints=cst, verbose=True)
    print(ss)
#print(s.eval())
#print(s.get_foldability())
