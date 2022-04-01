#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example::

    rna_mq_dfire.py AA/*.pdb
    id,fn,dfire
    1,1msy_A-000017_AA.pdb,-16181.328803
    2,1msy_A-000176_AA.pdb,-18172.496594
    3,1msy_A-000205_AA.pdb,-15526.172071

Output will be printed to stdout.
"""
from __future__ import print_function
import argparse
from rna_tools.tools.mq.Dfire.Dfire import Dfire
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument('-', "--", help="", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('-d', "--done", help="a csv with already done scores", default="")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    t = 'id,fn,dfire'
    print(t)
    t += '\n'
    for i, f in enumerate(args.file):
            wrapper = Dfire()
            result = wrapper.run(f, args.verbose)
            tl = ','.join([str(i + 1), os.path.basename(f), str(result)])
            print(tl)
            t += tl + '\n'

    with open('dfire.csv', 'w') as f:
        f.write(t)
