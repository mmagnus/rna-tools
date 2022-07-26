#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_calc_rmsd_all_vs_all.py - calculate RMSDs all vs all and save it to a matrix

Examples::

    rna_calc_rmsd_all_vs_all.py -i test_data -o test_output/rmsd_calc_dir.tsv
     # of models: 4
    ... 1 test_data/struc1.pdb
    ... 2 test_data/struc2.pdb
    ... 3 test_data/struc3.pdb
    ... 4 test_data/struc4.pdb

The program is using (https://github.com/charnley/rmsd).

You can also use PyMOL to do align or fit::

    python rna_calc_rmsd_all_vs_all.py -i test_data -o test_output/rmsd_calc_dir_align.mat -m align
     # of models: 5
    # test_data/2nd_triplex_FB_1AUA3_rpr.pdb test_data/struc1.pdb test_data/struc2.pdb test_data/struc3.pdb test_data/struc4.pdb
    0.0 4.13 4.922 4.358 4.368
    4.13 0.0 11.092 4.707 3.46
    4.922 11.092 0.0 11.609 11.785
    4.358 4.707 11.609 0.0 2.759
    4.368 3.46 11.785 2.759 0.0
    matrix was created!  test

"""
from __future__ import print_function

from rna_tools.tools.rna_calc_rmsd.lib.rmsd.calculate_rmsd import rmsd, get_coordinates, centroid, kabsch_rmsd
from rna_tools.tools.rna_calc_rmsd.rna_calc_rmsd import calc_rmsd_pymol

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
ic.configureOutput(prefix='> ')

import argparse
import glob
import re
import os


def get_rna_models_from_dir(directory):
    models = []
    if not os.path.exists(directory):
        raise Exception('Dir does not exist! ', directory)
    files = glob.glob(directory + "/*.pdb")
    files_sorted = sort_nicely(files)
    for f in files_sorted:
        models.append(f)
    return models


def sort_nicely(l):
    """Sort the given list in the way that humans expect.
    http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
    """
    def convert(text): return int(text) if text.isdigit() else text

    def alphanum_key(key): return [convert(c) for c in re.split('([0-9]+)', key)]
    l.sort(key=alphanum_key)
    return l


def calc_rmsd(a, b):
    """Calc rmsd."""
    atomsP, P = get_coordinates(a, None, None, 'pdb', True)
    atomsQ, Q = get_coordinates(b, None, None, 'pdb', True)

    # Calculate 'dumb' RMSD
    # gnormal_rmsd = rmsd(P, Q)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    # if False:
    #    V = rotate(P, Q)
    #    V += Qc
    #    write_coordinates(atomsP, V)
    #    quit()

    return kabsch_rmsd(P, Q)


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', "--input-dir",
                        default='',
                        help="input folder with structures")

    parser.add_argument('-o', "--matrix-fn",
                        default='matrix.txt',
                        help="ouput, matrix")

    parser.add_argument('-m', "--method",
                        default='all-atom',
                        help="all-atom, pymol: align, fit")


    # parser.add_argument("-s", "--save",
    #                    action="store_true", help="")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    input_dir = args.input_dir
    matrix_fn = args.matrix_fn

    models = get_rna_models_from_dir(input_dir)

    # print('number of models:', len(models))

    f = open(matrix_fn, 'w')
    t = '# '
    for r1 in models:
        # print r1,
        t += os.path.basename(r1) + ' '
    # print
    t += '\n'

    c = 1
    avg = 0
    avgs = []
    for r1 in models:
        for r2 in models:
            if args.method == 'align':
                rmsd_curr = calc_rmsd_pymol(r1, r2, 'align')[0] # (0.0, 0)
            elif args.method == 'fit':
                rmsd_curr = calc_rmsd_pymol(r1, r2, 'fit')[0] # (0.0, 0)                
            else:
                rmsd_curr = calc_rmsd(r1, r2)
            t += str(round(rmsd_curr, 3)) + ','
            avg += round(rmsd_curr, 3)
            avgs.append(round(rmsd_curr, 3))
        # print('...', c, r1)
        c += 1
        t = t[:-1] # remove ending ,
        t += '\n'

    f.write(t)
    f.close()

    def avgf(lst):
        return round(sum(lst) / len(lst), 2)
    ic.disable()
    ic(avg)
    ic(avgf(avgs))
    ic(c)
    ic(avgs)

    #print('score')
    print(avgf(avgs))
    #print(t.strip())  # matrix

    if 0:
        if True:
            print('matrix was created! ', matrix_fn)
        else:
            print('matrix NOT was created!')
