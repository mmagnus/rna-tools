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

The program is using (https://github.com/charnley/rmsd)
"""
from __future__ import print_function

from rna_tools.tools.rna_calc_rmsd.lib.rmsd.calculate_rmsd import rmsd, get_coordinates, centroid, kabsch_rmsd

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

    # parser.add_argument("-s", "--save",
    #                    action="store_true", help="")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    input_dir = args.input_dir
    matrix_fn = args.matrix_fn

    models = get_rna_models_from_dir(input_dir)

    print(' # of models:', len(models))

    f = open(matrix_fn, 'w')
    t = '# '
    for r1 in models:
        # print r1,
        t += str(r1) + ' '
    # print
    t += '\n'

    c = 1
    for r1 in models:
        for r2 in models:
            rmsd_curr = calc_rmsd(r1, r2)
            t += str(round(rmsd_curr, 3)) + ' '
        print('...', c, r1)
        c += 1
        t += '\n'

    f.write(t)
    f.close()

    print(t.strip())  # matrix

    if True:
        print('matrix was created! ', matrix_fn)
    else:
        print('matrix NOT was created!')
