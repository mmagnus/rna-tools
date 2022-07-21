#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
trx_score_onCol.py t3-3_rpr_fx.pdb --growth ~/Desktop/trx/mutget/t1-3/t1-3_full2.csv
"""

import pandas as pd
import argparse
import os

from rna_tools.tools.triplexibility import trx
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--fig", help="", default="fig.py") # nargs='+')
    parser.add_argument("ref", help="", default="") # nargs='+')
    parser.add_argument("--growth", help="", default="") # nargs='+')
    parser.add_argument("-e", "--edge", help="", default= "Triple_cWW_cHS") # nargs='+')
    parser.add_argument("-t", "--threshold", help="if t not none then if score > t is 1 else 0", type=int) # nargs='+')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    exec(open(args.fig).read())
    
    rows = [] # seq, scores
    for a, b in zip(aa.split(','), bb.split(',')):
         if a != 'x':
             for c in cc.split(','):

                 row = []
                 cols = []

                 cols.append('seq')
                 cols.append('seql')

                 f2 = args.ref + "_" + a + b + c + ".pdb"
                 f3 = f2.replace('.pdb', '_rpr.pdb')

                 csvf = a + b + c + ".csv"
                 seqm = a + b + c
                 row.append(seqm)
                 row.append(seqm.lower())
                 
                 df = pd.read_csv(csvf, sep=',', index_col=False)#, index=False)
                 score = float(df['rmsd'][0])
                 ic(score, args.threshold)
                 if args.threshold:
                     if score > args.threshold:
                         score = 1
                         ic(1)
                     else:
                         score = 0
                 ic(score)
                 cols.append('rmsd')
                 row.append(score)

                 dfexemplar = df[df.fn.str.contains("exemplar")]
                 ic(dfexemplar)
                 score_exem = float(dfexemplar['rmsd'].head(1))
                 if args.threshold:
                     if score_exem > args.threshold:
                         score_exem = 1
                     else:
                         score_exem = 0
                 ic(score_exem)
                 #except:
                 #    score_exem = 10 # there is not even any exemplar
                 row.append(score_exem)
                 cols.append('rmsd_exemplary')
                 
                 dfedge = df[df.fn.str.contains(args.edge)]
                 ic(dfedge)
                 try:
                     score_edge = float(dfedge['rmsd'].head(1))
                 except TypeError: # 
                     score_edge = 10

                 if args.threshold:
                     if score_edge > args.threshold:
                         score_edge = 1
                     else:
                         score_edge = 0
                 ic(score_edge)
                 row.append(score_edge)
                 cols.append('rmsd_edge')

                 scores = []
                 score_westhof = trx.wscore(seqm, 'cWW_cHS', type='exists')#clashes')
                 row.append(score_westhof)
                 cols.append('exists')
                 
                 score_westhof2 = trx.wscore(seqm, 'cWW_cHS', type='instances')
                 row.append(score_westhof2)
                 cols.append('instances')
                 
                 score_westhof3 = trx.wscore(seqm, 'cWW_cHS', type='contacts_bo')#clashes')
                 row.append(score_westhof3)
                 cols.append('contacts')

                 #score_westhof = trx.wscore(seqm, 'cWW_cHS')
                 #sscores_westhof.append(score_westhof)

                 if args.growth:
                     cols.append('growthby')
                     df = pd.read_csv(args.growth, sep=',', index_col=False)
                     growth = df[df['seql'] == seqm.lower()]['growthby'].values[0]
                     row.append(growth)

                 rows.append(row)                 
             
    print(cols)
    print(rows)
    
    df = pd.DataFrame(rows, columns=cols)
    print(df)
    df.to_csv('_scores_.csv', index=False)
