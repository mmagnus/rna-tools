#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The output is save for a csv file, <fn>.csv
"""
from __future__ import print_function

import argparse
import tempfile
import os
import time


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    try:
        import __main__
        __main__.pymol_argv = ['pymol', '-qc']
        import pymol  # import cmd, finish_launching
        pymol.finish_launching()
    except ImportError:
        print('calc_rmsd_pymol: you need to have installed PyMOL')
        sys.exit(0) # no error

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        t = tempfile.NamedTemporaryFile(delete=False) # True)

        pymol.cmd.reinitialize()
        pymol.cmd.delete('all')
        bn = f.split(os.sep)[-1]
        pymol.cmd.load(f, bn)
        pymol.cmd.do('contacts *')
        y = pymol.cmd.do("get_raw_distances contacts_polar, filename=" + t.name)

        time.sleep(1)
        nohb = len(open(t.name).readlines())

        with open(bn.replace('.pdb', '.csv'), 'w') as f:
           f.write(bn + ',' + str(nohb) + '\n')
        print(bn + ',' + str(nohb))
