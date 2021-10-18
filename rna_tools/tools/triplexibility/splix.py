#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
"""
def here which triple to chagne
"""

import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("seq", help="U57:A27.A53 / wt UAA", default="UAA") # nargs='+')
    parser.add_argument("--d2", help="U6-A59*U2-U23 / wt AU", default="UA") # nargs='+')
    parser.add_argument("--growth", help="experimental results", default=0) # nargs='+')
    return parser

def duplex_energy(s):
    """
    cGAUCgaaaGAUCg
    (((((....)))))
    rna_secondary_structure_prediction.py --method mcfold --file wt.fa
    """
    if s == 'AU':
        return -15.81
    if s == 'AA':
        return -13.42
    if s == 'UA':
        return -15.81
    if s == 'CA':
        return -13.0
    if s == 'GA':
        return -13.73
    # write them down
    
def duplex2_energy(s):
    s  = 'c' + s[1] + 'GCgaaaGC' + s[0] + 'c'
    ss = '((((....))))'
    from rna_tools.Seq import RNASequence
    seq = RNASequence(s)
    en, ss, comment = seq.predict_ss('mcfold', constraints=ss)
    print('   ', s, en, ss)    
    
def get_trx(seq, triple):
    # ic('>>', seq, triple)
    if triple == 't1_3':
        f = 'db/t1-3_cWW_tHW_UAA_exemplar_rpr_rmsd.csv'
    elif triple == 't2_3':
        f = 'db/t2-3-UAU_rmsd.csv'
    elif triple == 't2_4':
        f = 'db/t2-4-ACA_rmsd.csv'
    elif triple == 't3_3':
        f = 'db/t3-3_AUA_rmsd.csv'
    for l in open(f):
        y = -1
        if seq.lower() in l.lower():
            f, r = l.split(',') # Triple_cWW_tHS_AUG_exemplar_rpr.pdb,3.804
            r = float(r)
            #if r > 2: # 1.75:
            y = r
            print('   ' + l.strip())
            break
    return y, l.split(',')[0]

import pandas as pd
df = pd.read_csv('triple-db.csv')

def get_clashscore(t):
    """
    Triple_cWW_tHW_UAA_exemplar_rpr.pdb
    """
    t = t.replace('Triple_', '').replace('_rpr.pdb', '').replace('_exemplar', '')
    #ic(df)
    x = df[df['triple'] == t]
    ic(x)
    return int(x['clashes'])

def get_instances(t):
    """
    Triple_cWW_tHW_UAA_exemplar_rpr.pdb
    """
    t = t.replace('Triple_', '').replace('_rpr.pdb', '').replace('_exemplar', '')
    #ic(df)
    x = df[df['triple'] == t]
    # ic(x)
    return int(x['instances'])


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    s = args.seq
    d2 = args.d2
    
    from icecream import ic
    import sys
    ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))

    print('-' + s + '------------------------------------')
    sc = 0
    i = ic

    d = s[1] + s[0] # duplex a27:u57
    ic(d, duplex_energy(d))
    sc += 15.81 + duplex_energy(d) # 3 is a ad score

    ic(d)    
    duplex2_energy(d2[0] + d2[1])

    t1_3 = d2[0] + d2[1] + s[2] # U:A*A53
    ic(t1_3)#, get_trx(t1_3, 't1_3')))
    rmsd, triple = get_trx(t1_3, 't1_3')
    print(triple, rmsd)

    cs = get_clashscore(triple)
    inst = get_instances(triple)
    ic(cs, inst)

    t2_3 = d2[0] + d2[1] + s[0] # U2-A59:U2-U23*U57
    ic(t2_3) #, get_trx(t2_3, 't2_3'))
    rmsd, triple = get_trx(t2_3, 't2_3')
    #aaaa
    #msd, triple += get_trx(t2_3, 't2_3')
    cs = get_clashscore(triple)
    inst = get_instances(triple)
    ic(cs, inst)

    t2_4 = s[1] + 'C' + s[2] # ? # a27/c29/a53
    ic(t2_4)#, get_trx(t2_4, 't2_4'))
    rmsd, triple = get_trx(t2_4, 't2_4')
    cs = get_clashscore(triple)
    inst = get_instances(triple)
    ic(cs, inst)
    
    t3_3 = s[1] + s[0] + s[2] # a27/u57/53 ?
    ic(t3_3)#, get_trx(t3_3, 't3_3'))
    rmsd, triple = get_trx(t3_3, 't3_3')
    cs = get_clashscore(triple)
    inst = get_instances(triple)
    ic(cs, inst)

    g = float(args.growth) # ACG-0.5
    sc = round(sc, 2)
    if g == 1:
        print('^ grow:', g,  sc)
    elif g == 0:
        print('x dead ', g, sc)
    else:
        print('/ semi ', g, sc)
    # ic(t3_3, get_trx(t3_3, 't3_3'))
