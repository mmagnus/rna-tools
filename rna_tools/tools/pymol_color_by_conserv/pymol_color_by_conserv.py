#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
from Bio import SeqIO


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("conserv", help="consigns file", default="")
    parser.add_argument("alignment", help="", default="")
    parser.add_argument("chain", help="", default="A")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    consens = open(args.conserv).readline().strip().split(',')[1:]
    #print('Consensus:', consens)
    chain = args.chain

    for record in SeqIO.parse(args.alignment, "fasta"):
        header = record.id
        seq = record.seq
        break

    c = 1
    dat = []
    txt = ''
    seqq = ''
    bar = ''
    color = ''
    for s, cons in zip(seq, consens):
        if s != '-':
            dat.append([c, cons])
            seqq += s
            if float(cons) == 100:
                bar += '*'
                color = 'red'
            elif float(cons) > 90:
                color = 'red'
                bar += '▇'
            elif float(cons) > 75:
                color = 'orange'
                bar += '▆'
            elif float(cons) > 60:
                color = 'brown'
                bar += '▅'
            elif float(cons) > 45:
                color = 'yellow'
                bar += '▄'
            elif float(cons) > 30:
                color = 'gray30'
                bar += '▂'
            else:
                color = 'gray'
                bar += '▁'
            if not color:
              print(color, c, chain, cons)
            l = 'color %s, resi %s and chain %s; #%s' % (color, c, chain, cons)
            print(l)
            txt += l
            #print(c, s, cons)
            c += 1
    print('# ' + seqq)
    print('# ' + bar)
