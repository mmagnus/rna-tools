#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
proof of concept
"""
import pandas as pd

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
ic.configureOutput(prefix='> ')

import argparse
import glob
import random

TRX_DB = "/Users/magnus/work/src/rna-tools/rna_tools/tools/triplexibility/trx-db/triple-db-ext-2+Clash200x2+Contacts-1_fix36UCG.csv"


def get_seqs(type='all'):
    seq = []
    if type == 'all':    
        for a in ['a', 'c', 'g', 'u']:
            for b in ['a', 'c', 'g', 'u']:
                 for c in ['a', 'c', 'g', 'u']:
                     seq.append(a + b + c)
    if type == 't11':
        aa='x,G,G,c' # U2 G21
        bb='x,C,g,g' # U6 C61
        cc='U,c,a,g' # U6 U80
        # U2 G21 / U6 C61 / U6 U80
        for a, b in zip(aa.split(','), bb.split(',')):
             if a != 'x':
                 for c in cc.split(','):
                     seq.append(a + b + c)
        # C61g -> GgU
        seq = seq + ['GaU', 'GgU', 'GuU',
               'uaU', 'ugU', 'uuU',
               'aaU', 'agU', 'auU',
               'caU', 'cgU', 'cuU',       
               ]

    if type == 't12':
        aa='x,c,c,a,c'
        bb='x,g,u,u,c'
        cc='g,u,a,c'
        for a, b in zip(aa.split(','), bb.split(',')):
             if a != 'x':
                 for c in cc.split(','):
                     seq.append(a + b + c)

    return seq

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-d", "--debug",
                        action="store_true", help="be verbose")
    parser.add_argument("seq", help="")#, nargs='+')
    parser.add_argument("--edge", help="e.g., cWW_cHS")
    return parser


def growth_t13(s):#'uaa'):
    df = pd.read_csv('/Users/magnus/Desktop/trx/mutget/t1-3/t1-3_full_Yusuf_SimCols.csv')
    df = df.set_index('id')
    s = s.lower()
    # print(df)
    p = s[:2]
    t = s[2]
    return df[p][t]

def get_growth(s, t='t11'):
    """Add constrain on G for t11"""
    s = s.lower()
    if t == 't11':
        df = pd.read_csv('/Users/magnus/Desktop/trx/mutget/t1-1_scores_d8.csv')
        #ic(s)
        return int(df[df['seql'] == s]['growth'])

    if t == 't12':
        f = 'file:///Users/magnus/Desktop/trx/mutget/t1-2/t1-2/_scores_.csv'
        df = pd.read_csv(f)
        return int(df[df['seql'] == s]['growth'])
    
    if t == 't13':
        df = pd.read_csv('/Users/magnus/Desktop/trx/mutget/t1-3/t1_3_rmNan.csv')#t1_3.csv')
        return int(df[df['seql'] == s]['growthb'])
        
    df = df.set_index('id')
    p = s[:2]
    t = s[2]
    return df[p][t]

def growth_t33(s):#'uaa'):
    df = pd.read_csv('/Users/magnus/Desktop/trx/mutget/t3-3/fig_SimpCols.csv')
    df = df.set_index('id')
    s = s.lower()
    # print(df)
    p = s[:2]
    t = s[2]
    return df[p][t]

def growth_t24(s):
    df = pd.read_csv('/Users/magnus/Desktop/trx/mutget/t2-4/fig_SimpCols.csv')
    df = df.set_index('id')
    s = s.lower()
    # print(df)
    p = s[:2]
    t = s[2]
    return df[p][t]

def growth_t23(s):
    """Return -1 for Nan otherwise return value"""
    df = pd.read_csv('/Users/magnus/Desktop/trx/mutget/t2-3/t2-3_SimpCols.csv')
    df = df.set_index('id')
    s = s.lower()
    # print(df)
    p = s[:2]
    t = s[2]
    if pd.isna(df[p][t]):
        return -1
    else:
        return df[p][t]

def hr(t):
    l = len(t)
    half = int(l/2)
    print('-' * (40 - half) + ' ' + t + ' ' + '-' * (40 - half))

def fopen(filename, verbose=0):
    import subprocess
    cmd = 'mdfind -name "' + filename + '"'
    if verbose:
        print(cmd)
    out = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).stdout.read().decode()
    first_hit = out.split('\n')[0]
    return open(first_hit)

def tdb_score(seq, edge='cWW_cHW', type='instances', binary=False, if_exists=1, triplex=''):
    """
    x is::

        triple  instances  clashes  near  exists  contacts_bo  contacts_no_bk  contacts_rpr  clashscore_rpr  clashscore_no_bk  clashscore_bo  Unnamed: 11  Unnamed: 12  Unnamed: 13
        1557  cWW_cHS_UCU        1.0      0.0   0.0     1.0          2.0             2.0           3.0             0.0               0.0            0.0          NaN          NaN          NaN
        ('cWW_cHS_GCU:', 4) ('cWW_cHS_UCU:', 2)

        ValueError

               triple  instances  clashes  near  exists  contacts_bo  contacts_no_bk  contacts_rpr  clashscore_rpr  clashscore_no_bk  clashscore_bo  Unnamed: 11  Unnamed: 12  Unnamed: 13
    5202  cWW_cHS_CGU        0.0      0.0   0.0     0.0          NaN             NaN           NaN             NaN               NaN            NaN          NaN          NaN          NaN
    Traceback (most recent call last):
      File "splix_trx2.py", line 76, in <module>
        score2 = score_triple(triple2nd, edge)
      File "splix_trx2.py", line 53, in score_triple
        score = triple + ':', round(int(x['contacts_bo']))
      File "/Users/magnus/miniconda2/envs/py37/lib/python3.7/site-packages/pandas/core/series.py", line 128, in wrapper
        return converter(self.iloc[0])
    """
    if triplex == 't11':
        if seq.endswith('g'):
            ic('g correction used for ', seq)
            return 0

    import pandas as pd
    df = pd.read_csv('/Users/magnus/work/src/rna-tools/rna_tools/tools/triplexibility/triple-db.csv')
    triple = edge + '_' + seq.upper()
    #ic(triple)
    x = df[df['triple'] == triple]
    if if_exists:
        if not int(x['exists']):
            return -1
    # ic(seq, x)
    #r3 = fopen('Triple_' + edge + '_' + seq  + '_rpr.pdb.3dcnn.csv').read()
    #score = triple + ':', round(int(x[contacts_bo'instances']) / 18, 2)
    # print(x[type])
    try:
        score = str(round(int(x[type]))) # triple + ':' +  # contacts_bo
    except ValueError:
        score = -1 # triple + ':-1'
    #r3.split(',')[1].strip()
    score = int(score)
    #print(score)
    if binary:
        if score > 0:
            return 1
        else:
            return 0
    else:
        return score

def wall_score(seq, binary=False):
    """
    """
    import pandas as pd
    df = pd.read_csv("/Users/magnus/work/src/rna-tools/rna_tools/tools/triplexibility/db/triple-db-ext-2+Clash200x2+Contacts-1_fix36UCG.csv")
    seq = seq.upper()
    print(seq)
    x = df[df['seq'] == seq]
    #print(df)
    #ic(x)
    #r3 = fopen('Triple_' + edge + '_' + seq  + '_rpr.pdb.3dcnn.csv').read()
    #score = triple + ':', round(int(x[contacts_bo'instances']) / 18, 2)
    return x

def score_triple2(seq, edge, rand=False):
    """
    """
    x = df[df['seq'] == seq]
    #r3 = fopen('Triple_' + edge + '_' + seq  + '_rpr.pdb.3dcnn.csv').read()
    #score = triple + ':', round(int(x[contacts_bo'instances']) / 18, 2)
    try:
        score = str(round(int(x[edge + '_contacts_bo']))) # triple + ':' + 
    except:
        # print(x)
        score = -1 # triple + ':-1'
    #r3.split(',')[1].strip()
    score = int(score)
    binary = True
    if rand:
        return random.choice([0, 1])
    if binary:
        if score >= 0:
            return 1
        else:
            return 0

def _get_growth_(seq, dfg, a):
    """This is ugly, dfg global used"""
    x = dfg[dfg['seq'] == seq]
    try:
        g = int(x['growth' + a.upper()]) # , int(x['growthA27'])
    except:
        print(seq)
    binary = True
    if binary:
        if g > 0:
            return 1
        else:
            return 0

def init():
    import pandas as pd
    df = pd.read_csv("/Users/magnus/Desktop/trx/db/triple-db.csv")
    return df

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    edge = args.edge
    seq = args.seq
    v = args.verbose
    
    # if v: print(seq)
    # dfg = pd.read_csv("/Users/magnus/Desktop/trx/db/triple-db.csv")

    #edge1st = 'cWW_cHS'
    #triple1st = edge1st + '_' + 'GC' + seq[-1]
    #triple2nd = edge + '_' + seq

    #triple1st = seq
    # first
    #if v: print(triple1st)
    score = wscore(args.seq, args.edge)#, df)
    print(score)


    # second
    #if v: print(triple2nd)
    #score2 = score_triple(triple2nd, edge)
    #print(score, score2, get_growth(seq))
    sys.exit(1)

    x = growth_t23('AGU') # Nan
    ic(type(x))
    if x == 'nan':
        ic(x)
    import numpy as np
    if x == np.nan:
        print(x)
        ic(np.nan)

    wall_score('UAU')
