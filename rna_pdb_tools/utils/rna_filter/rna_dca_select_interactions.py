#!/usr/bin/env python
"""select top interactions"""

import pandas as pd
import sys


def sortdata(filename):
    df = pd.read_csv(filename, sep=" ", names=['i', ' ', 'j', ' ', ' ', 'scores'])
    df2 = df.sort_values(by=['scores'], ascending=False)
    maxvalue = max(df2[['j']].values)
    l = int(maxvalue)
    print 'l', l
    L = l / 2
    #L = l
    print 'L', L
    df3 = df2[0:L]
    df4 = df3[['i', 'j', 'scores']]
    df4['i'] -= 1
    df4['j'] -= 1
    print df4
    print 'Output file created:' + filename + "_selected.csv"
    df4.to_csv(filename + "_" + str(L) + "_sselected.csv", sep=" ")


def usage():
    if len(sys.argv) != 2:
        print sys.argv[0] + " <file> # is for sorting data"
        sys.exit()


if __name__ == "__main__":
    usage()
    sortdata(sys.argv[1])
