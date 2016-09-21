#!/usr/bin/env python

import pandas as pd
import sys

def sortdata(filename):
    df=pd.read_csv(filename,sep=" ",names=['i',' ','j',' ',' ','scores'])
    df2=df.sort_values(by=['scores'],ascending=False)
    maxvalue=max(df2[['j']].values)
    l=int(maxvalue)
    L=l/2
    print 'L', L
    df3=df2[0:L]
    df4=df3[['i', 'j', 'scores']]
    print df4
    print 'Output file created:' + filename + "_parsed.csv"
    df4.to_csv(filename+"_parsed.csv", sep=" ")

def usage():
    if len(sys.argv)!=2:
        print sys.argv[0] + " <file> # is for sorting data"
        sys.exit()

if __name__=="__main__":
    usage()
    sortdata(sys.argv[1])
    

 


