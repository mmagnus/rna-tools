import sys
from icecream import ic

ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')

import pandas as pd
df = pd.read_csv('triple-db.csv')
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")
 
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    #parser.add_argument("seq", help="U57:A27.A53 / wt UAA", default="UAA") # nargs='+')
    #parser.add_argument("--d2", help="U6-A59:U2-U23 / wt AU", default="UA") # nargs='+')
    #parser.add_argument("--growth", help="experimental results", default=0) # nargs='+')
    return parser

def score_trx_double(seq, tindex):
    rmsd, triple, rmsd2, triple2 = get_trx(seq, tindex) # t1_3, 't1_3')

    cs = get_clashscore(triple)
    inst = get_instances(triple)

    cs2 = get_clashscore(triple2)
    inst2 = get_instances(triple2)

    print('>> triplex 1 triple 3: \n   %s with rmsd score of %.2f of clashscore %i and instances %i' % (triple, rmsd, cs, inst))
    print('                        \n  %s with rmsd score of %.2f of clashscore %i and instances %i' % (triple2, rmsd2, cs2, inst2))
    return (triple, rmsd, cs, inst)

def score_trx(seq, triple_id):
    """
    seq: 'UUA'
    triple_id: 't1_3'
    """
    return get_trx(seq, triple_id)
    rmsd, triple, l = get_trx(seq, triple_id)
    ic(triple)
    cs = get_clashscore(triple)
    inst = get_instances(triple)

    print(triple_id + ': %s with rmsd score of %.2f of clashscore %i and instances %i' % (triple, rmsd, cs, inst))
    return (triple, rmsd, cs, inst)

def duplex_energy(s):
    """
    cGAUCgaaaGAUCg
    (((((....)))))
    rna_secondary_structure_prediction.py --method mcfold --file wt.fa
    """
    s  = 'cG' + s[1] + 'GCgaaaGC' + s[0] + 'Cg'
    ss = '(((((....)))))'
    from rna_tools.Seq import RNASequence
    seq = RNASequence(s)
    en, ss, comment = seq.predict_ss('mcfold', constraints=ss)
    print('   ', s, en, ss)
    return en
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
    return -1 # write them down
"""    
def duplex2_energy(s):
    s  = 'c' + s[1] + 'GCgaaaGC' + s[0] + 'c'
    ss = '((((....))))'
    from rna_tools.Seq import RNASequence
    seq = RNASequence(s)
    en, ss, comment = seq.predict_ss('mcfold', constraints=ss)
    # print('   ', s, en, ss)
    return en
    
def get_trx_double(seq, triple):
    # ic('>>', seq, triple)
    if triple == 't1_3':
        f = 'db/t1-3_cWW_tHW_UAA_exemplar_rpr_rmsd.csv'
    elif triple == 't2_3':
        f = 'db/t2-3-UAU_rmsd.csv'
    elif triple == 't2_4':
        f = 'db/t2-4-ACA_rmsd.csv'
    elif triple == 't3_3':
        f = 'db/t3-3_AUA_rmsd.csv'

    rmsd1, line1 = -1, -1
    rmsd2, line2 = -1, -1

    for l in open(f):
        if seq.lower() in l.lower():
            triple, rmsd = l.split(',') # Triple_cWW_tHS_AUG_exemplar_rpr.pdb,3.804
            if rmsd1 == -1:
                rmsd1 = round(float(rmsd), 2)
                line1 = triple
            else:
                rmsd2 = round(float(rmsd), 2)
                line2 = triple
                break
            #if r > 2: # 1.75:
            # print('   ' + l.strip())
    return rmsd1, line1, rmsd2, line2

def get_trx(seq, triple_id):
    """

    db/t1-3_cWW_tHW_UAA_exemplar_rpr_Triple_cWW_*_rpr.pdb_sugar_rmsd.csv'
       t1-3_cWW_tHW_UAA_exemplar_rpr_Triple_*ex*_rpr.pdb_sugar_rmsd.csv
    """
    
    rmsds = []
    lines = []
    scores = []
    headers = []

    c = 0 
    for s in ['t1-3_cWW_tHW_UAA_exemplar_rpr', 't2-3-UAU', 't2-4-ACA', 't3-3_AUA']:
        for mode in ['backbone+sugar']: # , 'c1', 'c1+Nx']: ['sugar', 
            for f in [#'Triple_cWW_*_rpr.pdb',
                      'Triple_*_rpr.pdb',
                      #'Triple_*ex*_rpr.pdb',
                      #'Triple_cWW_*ex*_rpr.pdb',
                      ]: #'Triple_*_rpr.pdb']:
                fdb = "db/%s_%s_%s_rmsd.csv" % (s, f.replace('*','X') ,mode)
                rmsd = -1
                if triple_id.replace('_', '-') in fdb:
                    for l in open(fdb):
                        if seq.lower() in l.lower() and rmsd == -1: # if -1 then there is no rmsd
                            # for this triple
                            triple, rmsd = l.split(',') # Triple_cWW_tHS_AUG_exemplar_rpr.pdb,3.804
                            cs = get_clashscore(triple)
                            if cs != 0:
                                rmsd = -1 # reset the flag
                                continue  # cs not zero, so go ahead
                            else:
                                ins = get_instances(triple)
                                rmsd = round(float(rmsd), 2)
                                headers.extend(['%s_%s_%s' % (s, f.replace('*','X') ,mode) + '_rmsd',
                                                '%s_%s_%s' % (s, f.replace('*','X') ,mode) + '_cs',
                                                '%s_%s_%s' % (s, f.replace('*','X') ,mode) + '_ins'])
                                scores.extend([rmsd, cs, ins])

                    if rmsd == -1: # still at the end of the loop
                            headers.extend(['%s_%s_%s' % (s, f.replace('*','X') ,mode) + '_rmsd',
                                            '%s_%s_%s' % (s, f.replace('*','X') ,mode) + '_cs',
                                            '%s_%s_%s' % (s, f.replace('*','X') ,mode) + '_ins'])
                            scores.extend([-1, -1, -1])

    return scores, headers

def get_trx_(seq, triple):
    # ic('>>', seq, triple)
    if triple == 't1_3':
        f = 'db/t1-3_cWW_tHW_UAA_exemplar_rpr_rmsd.csv'
    elif triple == 't2_3':
        f = 'db/t2-3-UAU_rmsd.csv'
    elif triple == 't2_4':
        f = 'db/t2-4-ACA_rmsd.csv'
    elif triple == 't3_3':
        f = 'db/t3-3_AUA_rmsd.csv'

    rmsd1, line1 = -1, -1
    rmsd2, line2 = -1, -1

    for l in open(f):
        if seq.lower() in l.lower():
            triple, rmsd = l.split(',') # Triple_cWW_tHS_AUG_exemplar_rpr.pdb,3.804
            rmsd = round(float(rmsd), 2)
            cs = get_clashscore(triple)
            ic(triple, rmsd, cs)
            if cs != 0:
                continue  # cs not zero, so go ahead
            else:
                break # stop this cycle and leave ;-)
    return rmsd, triple, l

def get_clashscore(t):
    """
    Triple_cWW_tHW_UAA_exemplar_rpr.pdb
    """
    t = t.replace('Triple_', '').replace('_rpr.pdb', '').replace('_exemplar', '')
    t = t.replace('tsS','tSS')
    x = df[df['triple'] == t]
    return int(x['clashes'])

def get_instances(t):
    """
    Triple_cWW_tHW_UAA_exemplar_rpr.pdb
    """
    t = t.replace('Triple_', '').replace('_rpr.pdb', '').replace('_exemplar', '')
    t = t.replace('tsS','tSS')
    x = df[df['triple'] == t]
    return int(x['instances'])


