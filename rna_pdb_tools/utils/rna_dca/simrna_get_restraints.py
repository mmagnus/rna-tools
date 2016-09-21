#!/usr/bin/python

"""Get the restrains for SimRNA based on a file with interactions.

simrna_get_restraints.py <interactions_mapped.txt>

Warning: it only works for chain A!
"""

import pandas as pd
import sys

def get_restraints(filename):
    restraints=[]
    df=pd.read_csv(filename,sep=" ")
    txt = ''
    for row in df.iterrows():
        i = row[1]["i"]
        j = row[1]["j"]
        txt += 'SLOPE A/' + str(i) + '/MB A/' + str(j) + '/MB 3 8 0.5 \n' # energy values about restrains. slope <-> penalty
        txt += 'WELL A/' + str(i) + '/MB A/' + str(j) + '/MB 3 8 0.5 \n' # well <-> reward # check paper
    print txt.strip()

def usage():
    if len(sys.argv)<2:
        print sys.argv[0] + " <interactions_parsed_mapped.csv> "
        print '           output goes to stdout'
        sys.exit()

if __name__=="__main__":
    usage()
    get_restraints(sys.argv[1])
