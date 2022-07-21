#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
"""
def here which triple to chagne
"""

from lib import *

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    from icecream import ic
    import sys
    ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
    ic.configureOutput(prefix='> ')

    sc = 0

    d2 = "UA"

    rows = []

    # encode seq
    
    
    #for t in ['t1_3', 't2_3', 't2_4', 't3_3']:
    #    headers.append(t + '_triple')
    #    headers.append(t + '_rmsd')
    #    headers.append(t + '_clashscore')
    #    headers.append(t + '_instances')        
    #headers.append('growth')

    for (s, growth) in [('UAA', 1), ('UGA', 0), ('UUA', 0), ('UCA', 0.75),
                        ('UAC', 1), ('UGC', 0), ('UUC', 0), ('UCC', 1),

                        ('AAC', 1), ('AGC', 1), ('AUC', 1), ('ACC', 1),
                        ('CAC', 1), ('CGC', 1), ('CUC', 0), ('CCC', 0.25),

                        ('UAG', 1), ('UGG', 0), ('UUG', 0), ('UCG', 0),
                        ('AAG', 0), ('AGG', 0), ('AUG', 0), ('ACG', 0.5),

                        ('CAG', 0), ('CGG', 0.75), ('CUG', 0), ('CCG', 0),
                        ('GAG', 0), ('GGG', 0.5), ('GUG', 0), ('GCG', 0.2),

                        ('UAU', 1), ('UGU', 0), ('UUU', 0), ('UCU', 0),
                        ('AAU', 1), ('AGU', 0.75), ('AUU', 0.25), ('ACU', 0.25),

                        ('CAU', 0), ('CGU', 0), ('CUU', 0), ('CCU', 0),
                        ('GAU', 0.25), ('GGU', 0.25), ('GUU', 0), ('GCU', 0),

                        ('AAA', 1), ('AGA', 1), ('AUA', 1), ('ACA', 1),
                        ('CAA', 1), ('CGA', 1), ('CUA', 1), ('CCA', 0.5),
                        ('GAA', 1), ('GGA', 1), ('GUA', 1), ('GCA', 1)]:
        print(s + ': ===========================================')
        d = s[1] + s[0] # duplex a27:u57

        row = [s[0], s[1], s[2]]
        headers = ['U57', 'A27', 'A53', 'd1_energy', 'd2_energy']  # 'seq', #['A', 'C', 'U', 'G',  

        if True:
            energy = duplex_energy(d)
            print('    calc energy of U2-A27:U6-U57 [AU] duplex now:', s[1] + s[0], energy)
            row.append(energy)
            #sc += 15.81 + duplex_energy(d) # 3 is a ad score

            energy = duplex2_energy(d2[0] + d2[1])
            print('    calc energy of U6-A59:U2-U23 [UA] duplex now:', d2[0] + d2[1], energy)
            d2_energy = energy
            row.append(energy)

        t1_3 = d2[0] + d2[1] + s[2] # U:A*A53
        t2_3 = d2[0] + d2[1] + s[0] # U2-A59:U2-U23*U57
        t2_4 = s[1] + 'C' + s[2] # ? # a27/c29/a53
        t3_3 = s[1] + s[0] + s[2] # a27/u57/53 ?

        headers += ['t1_3_1', 't1_3_2', 't1_3_3',
               't2_2_1', 't2_2_2', 't2_2_3',
               't2_4_1', 't2_4_2', 't2_4_3',
               't3_3_1', 't3_3_2', 't3_3_3']
               
        row += [t1_3[0], t1_3[1], t1_3[2],
               t2_3[0], t2_3[1], t2_3[2],
               t2_4[0], t2_4[1], t2_4[2],
               t3_3[0], t3_3[1], t3_3[2]]
               
        for triple_id, triple_id_t in [(t1_3, 't1_3'),
                                       (t2_3, 't2_3'),
                                       (t2_4, 't2_4'),
                                       (t3_3, 't3_3'),
                                       ]: 
            rowt, headerst = score_trx(triple_id, triple_id_t) # scores, headers
            row += rowt
            headers += headerst # header once is fine # ugly

        if 0:
            ic(len(row))
            ic(len(headers))

            c = 0
            for r, h in zip(row, headers):
                c += 1
                print(c, h, r)

        g = float(growth) # ACG-0.5

        ## sc = round(sc, 2)
        ## if g == 1:
        ##     print('^ grow:', g,  sc)
        ## elif g == 0:
        ##     print('x dead ', g, sc)
        ## else:
        ##     print('/ semi ', g, sc)

        # growth binary
        if growth > 0:
            row += [1]
        else:
            row += [0]
        # collects rows
        rows.append(row)

    headers.append('growth')

    with open('splix.csv', 'w') as f:
        f.write(','.join([str(h) for h in headers]) + '\n')
        for row in rows:
            f.write(','.join([str(x) for x in row]) + '\n')
        # ic(t3_3, get_trx(t3_3, 't3_3'))
        
