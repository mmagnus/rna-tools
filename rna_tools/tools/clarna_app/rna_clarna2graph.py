#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Get the input::

    rna_pdb_tools.py --extract 'A:80+A:61+B:21+A:85' S*.pdb --here
    for i in *extr.pdb; do echo $i; rna_clarna_run.py -ipdb $i; done | tee log.txt

Example::

    $ rna_clarna2graph.py cl.txt
      source target    n
    0     21     61  408
    1     80     85   24
    2     61     80    9
    > out: 'cl_short.csv'

and with ```--full```::

       source target     type    n
    0      21     61   WW_cis  405
    1      80     85   WW_cis   13
    2      80     85  HW_tran    6
    3      80     85  WH_tran    5
    4      61     80  WS_tran    4
    5      61     80   HS_cis    3
    6      21     61   WH_cis    1
    7      61     80   WS_cis    1
    8      21     61   SW_cis    1
    9      61     80   HW_cis    1
    10     21     61  SH_tran    1
    > out: 'cl_full.csv'

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-f', "--full",  action="store_true", help="split interactions into types")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    full = args.full
    v = args.verbose
    

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        print(f)
        all = {}
        c = 0
        for line in open(f):
            if line.strip():
                if v: ic(l.strip())
                l = line.split()
                if len(l) == 9:
                    # ['B', '21', 'A', '61', 'bp', 'U', 'A', 'WW_cis', '0.8365']
                    if full: # ('21-61-WW_cis', 405)
                        # k = l[1] + '-' + l[3] + '-' + l[7] # no resin
                        k = l[5] + l[1] + '-' + l[6] + l[3] + '-' + l[7]
                    else:
                        # k = l[1] + '-' + l[3] # only nodes ('21-61', 408)
                        k = l[5] + l[1] + '-' + l[6] + l[3]
                    if k in all:  # count endge
                        all[k] += 1
                    else:
                        all[k] = 1
            if '.pdb' in line.lower():
                c += 1
        # sort it
        import operator
        alls = sorted(all.items(), key=operator.itemgetter(1), reverse=True)
        nlst = []
        for i in alls:
            if v: ic(i)
            ne = i[0].split('-') + [i[1]]  # [['21', '61', 408], ['80', '85', 24], ['61', '80', 9]]
            nlst.append(ne)
        if v: ic(nlst)

        import pandas as pd
        if full:
            columns=['source', 'target', 'type', 'n']
        else:
            columns=['source', 'target', 'n']
        df = pd.DataFrame(nlst, columns = columns)
        df['n'] = df['n'] / c
        print('# %i' % c)
        print(df)
        if full:
            out = f.replace('.txt', '_full.csv')
        else:
            out = f.replace('.txt', '_short.csv')

        ic(out)
        df.to_csv(out, index=False)
