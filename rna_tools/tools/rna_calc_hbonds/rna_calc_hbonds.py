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
        from pymol import stored
        pymol.finish_launching()
    except ImportError:
        print('you need to have installed PyMOL')
        import sys
        sys.exit(0) # no error

    if list != type(args.file):
        args.file = [args.file]

    pymol.cmd.reinitialize()
    for index, f in enumerate(args.file):
        t = tempfile.NamedTemporaryFile(delete=False) # True)

        pymol.cmd.delete('all')
        bn = f.split(os.sep)[-1]
        pymol.cmd.load(f, bn)
        pymol.cmd.do('contacts *,*')
        pymol.cmd.save('out.pse')
        #pymol.cmd.do("print get_raw_distances('contacts_polar', filename='')")# , filename='" + t.name + "')")
        pymol.cmd.do("stored.y = get_raw_distances('contacts_polar_ok', filename='out.csv')",)# , filename='" + t.name + "')")
        # Triple_tWW_tSW_UAU_rpr.pdb,([(('Triple_tWW_tSW_UAU_rpr.pdb', 14), ('Triple_tWW_tSW_UAU_rpr.pdb', 31), 2.68421207936708), (('Triple_tWW_tSW_UAU_rpr.pdb', 37), ('Triple_tWW_tSW_UAU_rpr.pdb', 13), 2.714407209563896), (('Triple_tWW_tSW_UAU_rpr.pdb', 56), ('Triple_tWW_tSW_UAU_rpr.pdb', 33), 2.951731552332452)], [['Triple-U1948', 'Triple-A1915'], ['Triple-A1915', 'Triple-U1948'], ['Triple-U779', 'Triple-A1915']])
        nohb = len(stored.y[1]) # get to the list from the comment above
        # an old mechanism to open output file to see how many lines are there
        # nohb = len(open(t.name).readlines())
        with open(bn.replace('.pdb', '.csv'), 'w') as f:
            f.write(bn + ',' + str(nohb) + '\n')
        print(str(index + 1) + '\t' + str(round((index/len(args.file)) * 100)) + '%\t' + bn + '\t' + str(nohb))
