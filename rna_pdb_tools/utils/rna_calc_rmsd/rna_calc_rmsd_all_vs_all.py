#!/usr/bin/python
#-*- coding: utf-8 -*-
"""
rna_calc_rmsd git:(master) ✗ ./rna_calc_rmsd_all_vs_all.py
calc_rmsd_dir
--------------------------------------------------------------------------------
Usage: rna_calc_rmsd_all_vs_all.py [<options>]

Options:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir=INPUT_DIR
  -o MATRIX_FN, --matrix_fn=MATRIX_FN
                        ouput, matrix
  -s, --save

required biopython
"""
import Bio.PDB.PDBParser
import Bio.PDB.Superimposer
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO, Superimposer

from lib.rmsd.calculate_rmsd import *

import optparse
import sys
import math
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

def sort_nicely( l ):
   """ Sort the given list in the way that humans expect.

   http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
   """
   convert = lambda text: int(text) if text.isdigit() else text
   alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
   l.sort( key=alphanum_key )
   return l

def calc_rmsd(a,b):
    """empty selection"""
    atomsP, P = get_coordinates(a, None, None, 'pdb', True)
    atomsQ, Q = get_coordinates(b, None, None, 'pdb', True)

    # Calculate 'dumb' RMSD
    normal_rmsd = rmsd(P, Q)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    if False:
        V = rotate(P, Q)
        V += Qc
        write_coordinates(atomsP, V)
        quit()

    return kabsch_rmsd(P, Q)    

if __name__ == '__main__':
    print 'calc_rmsd_dir'
    print '-' * 80
    
    optparser=optparse.OptionParser(usage="%prog [<options>]")

    optparser.add_option('-i',"--input_dir", type="string",
                         dest="input_dir",
                         default='',
                         help="")

    optparser.add_option('-o',"--matrix_fn", type="string",
                         dest="matrix_fn",
                         default='matrix.txt',
                         help="ouput, matrix")

    optparser.add_option("-s", "--save",
                     action="store_true", default=False, dest="save", help="")

    
    (opts, args)=optparser.parse_args()

    if len(sys.argv) == 1:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)

    input_dir = opts.input_dir
    matrix_fn = opts.matrix_fn

    models = get_rna_models_from_dir(input_dir)        

    print ' # of models:', len(models)

    f = open(matrix_fn, 'w')
    t = '# '
    for r1 in models:
        #print r1,
        t += str(r1) + ' '
    #print
    t += '\n'

    c = 1
    for r1 in models:
            for r2 in models:
                rmsd_curr = calc_rmsd(r1, r2)
                t += str(round(rmsd_curr,3)) + ' '
            print '...', c, r1
            c += 1
            t += '\n'
            
    f.write(t)
    f.close()

    print t.strip() # matrix

    if True:
        print 'matrix was created! ', matrix_fn
    else:
        print 'matrix NOT was created!'
