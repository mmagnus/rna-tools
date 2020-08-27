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
    parser.add_argument("molecule", help="", default="")
    parser.add_argument("alignment", help="", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    consens = open(args.conserv).readline().strip().split(',')[1:]

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
    print('Name,Conservation')
    for s, cons in zip(seq, consens):
        if s != '-':
            con = float(cons.strip()) / 100
            print(args.molecule + '-' + s + str(c) + ',' + str(round(con, 2))) # CWC15-Q4
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
            c += 1
    print('# ' + seqq)
    print('# ' + bar)
