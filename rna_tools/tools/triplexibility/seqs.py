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
    parser.add_argument("seq", help="", default="U57:A27.A53/ UAA") # nargs='+')
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
    
def score(seq, triple):
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
    return y

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    s = args.seq

    from icecream import ic
    import sys
    ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))

    print('-' + s + '------------------------------------')
    sc = 0
    i = ic

    d = s[1] + s[0] # duplex a27:u57
    ic(d, duplex_energy(d))
    sc += 15.81 + duplex_energy(d) # 3 is a ad score

    t1_3 = 'UA' + s[2]
    ic(t1_3)#, score(t1_3, 't1_3')))
    sc += score(t1_3, 't1_3')

    t2_3 = 'UA' + s[0]
    ic(t2_3) #, score(t2_3, 't2_3'))
    sc += score(t2_3, 't2_3')

    t2_4 = 'AC' + s[2] # ?
    ic(t2_4)#, score(t2_4, 't2_4'))
    sc += score(t2_4, 't2_4')
    
    t3_3 = s[1] + s[0] + s[2] # a27/u57/53 ?
    ic(t3_3)#, score(t3_3, 't3_3'))
    sc += score(t3_3, 't3_3')

    g = float(s.split('-')[1]) # ACG-0.5
    sc = round(sc, 2)
    if g == 1:
        print('^ grow:', g,  sc)
    elif g == 0:
        print('x dead ', g, sc)
    else:
        print('/ semi ', g, sc)
    # ic(t3_3, score(t3_3, 't3_3'))
    
