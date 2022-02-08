#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--fig", help="", default="fig.py") # nargs='+')
    parser.add_argument("ref", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    exec(open(args.fig).read())
    
    sl = []
    sle = []
    sledge = []
    for a, b in zip(aa.split(','), bb.split(',')):
         if a != 'x':
             scores = []
             scores_exem = []
             scores_edge = []
             
             for c in cc.split(','):

                 f2 = args.ref + "_" + a + b + c + ".pdb"
                 f3 = f2.replace('.pdb', '_rpr.pdb')

                 csvf = a + b + c + ".csv"
                 seqm = a + b + c
                 
                 df = pd.read_csv(csvf, sep=',', index_col=False)#, index=False)
                 score = df['rmsd'][0]
                 print(score)
                 scores.append(score)

                 dfexemplar = df[df.fn.str.contains("exemplar")]
                 print(dfexemplar)
                 score_exem = float(dfexemplar['rmsd'].head(1))
                 print(score_exem)
                 #except:
                 #    score_exem = 10 # there is not even any exemplar
                 scores_exem.append(score_exem)
                 
                 dfedge = df[df.fn.str.contains("Triple_cWW_cHS")]
                 print(dfedge)
                 try:
                     score_edge = float(dfedge['rmsd'].head(1))
                 except TypeError: # empty
                     score_edge = 10
                 print(score_edge)
                 scores_edge.append(score_edge)

             sl.append(scores)
             sle.append(scores_exem)
             sledge.append(scores_edge)

    df = pd.DataFrame(sle)
    df = df.T
    df.columns = ['a', 'b', 'c']
    df.to_csv('sle.csv', index=False)

    df = pd.DataFrame(sledge)
    df = df.T
    df.columns = ['a', 'b', 'c']
    df.to_csv('sledge.csv', index=False)
