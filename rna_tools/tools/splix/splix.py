#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
import glob


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-d", "--debug",
                        action="store_true", help="be verbose")
    parser.add_argument("mutant", help="", nargs='+')
    return parser


def hr(t):
    l = len(t)
    half = int(l/2)
    print('-' * (40 - half) + ' ' + t + ' ' + '-' * (40 - half))
    
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    hr('Literature:')
    lit = open('literature.tsv').read()
    g = []
    for i in lit.split('\n'):
        if i.strip():
            x, y = i.split('	') # tab
            if x == y:
                continue
            e = [x, y]
            g.append(e)

    t = []
    for m in args.mutant:
        for n in g: # graph
            if m.lower() in n[1].lower():  #
            # if n[0].upper() == m.upper():
                #print(n[0].ljust(5) + ' ' + n[1])
                n1 = n[1].strip()
                if n1 not in t:
                    t.append(n1)
    for l in t:
        print(l)
        print()

    hr('Secondary structure:')
    g2 = [
        ["U2-A25", "U2-U10"],
       ]

    g = g2
    if args.debug: print(g)
    for m in args.mutant:
        for n in g:
            if n[0].upper() == m.upper():
                print(n[0].ljust(5) + ' ' + n[1])

    for f in glob.glob('db/*'):
        #g3 = """
        #U2-A25,CEF1-H31
        #U2-A25,PRP8-D1094
        #"""
        hr(f)
        g3 = open(f).read()
        s = []
        for i in g3.split('\n'):
            if i.strip():
                x, y = i.split(',')
                if x == y:
                    continue
                e = [x, y]
                if e not in s:
                    s.append(e)
                e = [y, x]
                if e not in s:
                    s.append(e)

        if args.debug: print(g3)
        for m in args.mutant:
            for n in s:
                if n[0].upper() == m.upper():
                    print(n[0].ljust(5) + ' ' + n[1])
            

    
