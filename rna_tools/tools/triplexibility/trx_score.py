#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
    parser.add_argument("-e", "--edge", help="", default= "Triple_cWW_cHS") # nargs='+')
    parser.add_argument("-t", "--threshold", help="if t not none then if score > t is 1 else 0", type=int) # nargs='+')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    exec(open(args.fig).read())
    
    sl = []
    sle = []
    sledge = []
    stdb = []
    stdb2 = []
    stdb3 = []

    for a, b in zip(aa.split(','), bb.split(',')):
         if a != 'x':
             scores = []
             scores_exem = []
             scores_edge = []
             scores_westhof2 = []
             scores_westhof3 = []
             scores_westhof = []
             
             for c in cc.split(','):

                 f2 = args.ref + "_" + a + b + c + ".pdb"
                 f3 = f2.replace('.pdb', '_rpr.pdb')

                 csvf = a + b + c + ".csv"
                 seqm = a + b + c
                 
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
                 scores.append(score)

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
                 scores_exem.append(score_exem)
                 
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
                 scores_edge.append(score_edge)

                 scores = []
                 score_westhof = trx.wscore(seqm, 'cWW_cHS', type='exists')#clashes')
                 scores_westhof.append(score_westhof)
                 score_westhof2 = trx.wscore(seqm, 'cWW_cHS', type='instances')
                 scores_westhof2.append(score_westhof2)
                 score_westhof3 = trx.wscore(seqm, 'cWW_cHS', type='contacts_bo')#clashes')
                 scores_westhof3.append(score_westhof3)

                 #score_westhof = trx.wscore(seqm, 'cWW_cHS')
                 #sscores_westhof.append(score_westhof)

             sl.append(scores)
             sle.append(scores_exem)
             sledge.append(scores_edge)

             stdb.append(scores_westhof)
             stdb2.append(scores_westhof2)
             stdb3.append(scores_westhof3)
             n = sl

    df = pd.DataFrame(sl)
    df = df.T
    df.to_csv('sl.csv', index=False)

    df = pd.DataFrame(sle)
    df = df.T
    df.to_csv('slexemplary.csv', index=False)

    df = pd.DataFrame(sledge)
    df = df.T
    df.to_csv('sledge.csv', index=False)

    df = pd.DataFrame(stdb)
    df = df.T
    df.to_csv('stdb.csv', index=False)

    df = pd.DataFrame(stdb2)
    df = df.T
    df.to_csv('stdb2.csv', index=False)

    df = pd.DataFrame(stdb3)
    df = df.T
    df.to_csv('stdb3.csv', index=False)

    ic(sle, stdb, stdb2, stdb3)
