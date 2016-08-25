#!/usr/bin/env python

import pandas as pd
import sys
import numpy as np

def rna_dca_mapping(fileseq,file_interactions):
    f = open(fileseq)
    header = f.readline().strip() # get rid of header
    seq = f.read().strip()
    print header
    print seq

    df=pd.read_csv(file_interactions,sep=" ")
    interactions = zip(df['i'].tolist(), df['j'].tolist())
    fileseq.strip() 
    print 'input:', fileseq.strip()
    print 'interactions:\n', interactions
    mapped_interactions = []
    v = True # verbose
    for i in interactions:
        ij = [0,0]
        ij[0] = i[0] - 1
        ij[1] = i[1] - 1
        if v: print seq[ij[0]], seq[ij[1]]
        if seq[ij[0]] == '.' or seq[ij[1]] == '.':
            if v: print 'Removed Interaction:', ij
        else:
            [a,b]=[ij[0] - seq[:ij[0]].count('.')+1, ij[1] - seq[:ij[1]].count('.')+1,]
            print i, '->', a,b
            mapped_interactions.append([a,b])

    print 'Mapped Interactions:\n', mapped_interactions
    print 'get_dists(' + str(mapped_interactions) + ')'
    for i in mapped_interactions:
        if v:print str(i)
    print 'output seq:\n', seq.replace('.','')
    print 'output file:', file_interactions+"_mapped.csv"

    a=pd.DataFrame(mapped_interactions, columns=["i","j"])
    a.to_csv(file_interactions+"_mapped.csv",sep=" ")

def usage():
    if len(sys.argv)!=3:
        print sys.argv[0] + " <seq.fa> <interactions_parsed.txt>"
        sys.exit()
if __name__=="__main__":
    usage()
    rna_dca_mapping(sys.argv[1],sys.argv[2])
