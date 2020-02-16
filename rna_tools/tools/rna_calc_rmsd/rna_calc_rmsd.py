#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""rna_calc_rmsd.py - calculate RMSDs of structures to the target

If you still have problem with various number of atoms, check out this issue: get_rnapuzzle_ready: how to deal with `Alternate location indicator (https://github.com/mmagnus/rna-pdb-tools/issues/30).

The program is using (https://github.com/charnley/rmsd).

Example #1::

    $ rna_calc_rmsd.py -t 6_0_solution_4GXY_rpr.pdb --model-selection=A:1-17+24-110+115-168 *.pdb
    rmsd_calc_rmsd_to_target
    --------------------------------------------------------------------------------
    method: all-atom-built-in
    # of models: 35
    6_0_solution_4GXY_rpr.pdb 0.0 3409
    6_Blanchet_1_rpr.pdb 22.31 3409
    6_Blanchet_2_rpr.pdb 21.76 3409
    6_Blanchet_3_rpr.pdb 21.32 3409
    6_Blanchet_4_rpr.pdb 22.22 3409
    6_Blanchet_5_rpr.pdb 24.17 3409
    6_Blanchet_6_rpr.pdb 23.28 3409
    6_Blanchet_7_rpr.pdb 22.26 3409
    6_Bujnicki_1_rpr.pdb 36.95 3409
    6_Bujnicki_2_rpr.pdb 30.9 3409
    6_Bujnicki_3_rpr.pdb 32.1 3409
    6_Bujnicki_4_rpr.pdb 32.04 3409
    ...

Example #2::

    time rmsd_calc_to_target.py
      -t 5k7c_clean_onechain_renumber_as_puzzle_srr.pdb
      --target-selection A:1-48+52-63
      --model-selection A:1-48+52-63
      --target-ignore-selection A/57/O2\\'
      clusters/*_AA.pdb

    rmsd_calc_rmsd_to_target
    --------------------------------------------------------------------------------
      target_selection:  A:1-48+52-63
      model_selection:   A:1-48+52-63
      target_ignore_selection:  A/57/O2'
      model_ignore_selection:
    # of models: 801
    fn,rmsd_all
    pistol_thrs0.50A_clust01-000001_AA.pdb,7.596
    pistol_thrs0.50A_clust02-000001_AA.pdb,7.766
    pistol_thrs0.50A_clust03-000001_AA.pdb,18.171
    [..]
    pistol_thrs0.50A_clust799-000001_AA.pdb,5.356
    pistol_thrs0.50A_clust800-000001_AA.pdb,15.282
    pistol_thrs0.50A_clust801-000001_AA.pdb,16.339
    # of atoms used: 1237
    csv was created!  rmsds.csv
    rmsd_calc_to_target.py -t 5k7c_clean_onechain_renumber_as_puzzle_srr.pdb
    37.93s user 1.07s system 87% cpu 44.650 total

"""
from __future__ import print_function

from rna_tools.tools.rna_calc_rmsd.lib.rmsd.calculate_rmsd import *
import sys
from rna_tools.tools.extra_functions.select_fragment import select_pdb_fragment_pymol_style, select_pdb_fragment
import argparse
import sys
import math
import glob
import re
import os
import pandas as pd


def get_rna_models_from_dir(files):
    """
    :param models: a list of filenames

    Example of the list::

       ['test_data/rp17/2_restr1_Michal1.pdb_clean.pdb', 'test_data/rp17/2a_nonrestr2_Michal1.pdb_clean.pdb',
       'test_data/rp17/3_nonrestr1_Michal1.pdb_clean.pdb', 'test_data/rp17/5_restr1_Michal3.pdb_clean.pdb']"""

    models = []
    #if not os.path.exists(directory):
    #    raise Exception('Dir does not exist! ', directory)
    #files = glob.glob(directory + "/*.pdb")
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

def calc_rmsd_pymol(pdb1, pdb2, method):
    """Calculate rmsd using PyMOL. Two methods are available: align and fit

    See:

    -  Align: http://www.pymolwiki.org/index.php/Align
    -  Fit:   http://www.pymolwiki.org/index.php/Fit

    Align can return a list with 7 items:

    RMSD after refinement
    Number of aligned atoms after refinement
    Number of refinement cycles
    RMSD before refinement
    Number of aligned atoms before refinement
    Raw alignment score
    Number of residues aligned

    in this version of function, the function returns `RMSD before refinement`.

    Install on OSX: ``brew install homebrew/science/pymol`` and set ``PYTHONPATH`` to
    your PyMOL packages, .e.g ::

      PYTHONPATH=$PYTHONPATH:/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages

    If problem::

      Match-Error: unable to open matrix file '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/data/pymol/matrices/BLOSUM62'.

    then define ``PYMOL_PATH`` in your .bashrc, e.g.::

       export PYMOL_PATH=/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/pymol/
     """

    try:
        import __main__
        __main__.pymol_argv = ['pymol', '-qc']
        import pymol  # import cmd, finish_launching
        pymol.finish_launching()
    except ImportError:
        print('calc_rmsd_pymol: you need to have installed PyMOL')
        sys.exit(0) # no error

    pymol.cmd.reinitialize()
    pymol.cmd.delete('all')
    pymol.cmd.load(pdb1, 's1')
    pymol.cmd.load(pdb2, 's2')
    if method == 'align':
        # experiments with align <https://pymolwiki.org/index.php/Align>
        # quiet = 0/1: suppress output {default: 0 in command mode, 1 in API}
        return  (pymol.cmd.align('s1', 's2',quiet=1, object='aln')[3],0) #, pymol.cmd.align('s1','s2')[4])
        #raw_aln = pymol.cmd.get_raw_alignment('aln')
        #print raw_aln
        #for idx1, idx2 in raw_aln:
        #    print '%s`%d -> %s`%d' % tuple(idx1 + idx2)
        #pymol.cmd.save('aln.aln', 'aln')

    if method == 'fit':
        return (pymol.cmd.fit('s1', 's2'), 'fit')

def calc_rmsd(a,b, target_selection, target_ignore_selection, model_selection, model_ignore_selection, verbose):
    """
    a is model
    b is target

    :params: a = filename of structure a
    :params: b = filename of structure b

    :return: rmsd, number of atoms
    """
    if verbose: print('in:', a)
    atomsP, P = get_coordinates(a, model_selection, model_ignore_selection, 'pdb', True)
    atomsQ, Q = get_coordinates(b, target_selection,target_ignore_selection,  'pdb', True)

    if verbose:
        print(atomsP, P)
        print(atomsQ, Q)

    if atomsQ != atomsP:
        print('Error: # of atoms is not equal target (' + b + '):' + str(atomsQ) + ' vs model (' + a + '):' + str(atomsP))
        return (-1,0) # skip this RNA
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

    return round(kabsch_rmsd(P, Q),2), atomsP

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)#formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t',"--target-fn",
                            default='',
                         help="pdb file")

    parser.add_argument("--target-selection",
                            default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")

    parser.add_argument("--target-ignore-selection",
                            default='',
                         help="A/10/O2\'")

    parser.add_argument("--model-selection",
                            default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")

    parser.add_argument("--model-ignore-selection",
                            default='',
                         help="A/10/O2\'")

    parser.add_argument('-m', "--method",
                         default='all-atom-built-in',
                         help="align, fit")

    parser.add_argument('-o', "--rmsds-fn",
                         default='rmsds.csv',
                         help="ouput, matrix")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose")

    parser.add_argument('-pr', '--print-results',
                         action="store_true")

    parser.add_argument('-sr', '--sort-results',
                         action="store_true")

    parser.add_argument('-pp', '--print-progress',
                         default=False,
                         action="store_true")

    parser.add_argument("--target-column-name", action="store_true",
                        help="")

    parser.add_argument('files', help='files', nargs='+')

    return parser
# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    input_files = args.files  # opts.input_dir
    rmsds_fn = args.rmsds_fn
    target_fn = args.target_fn
    method = args.method

    print('# method:', method)

    target_selection = select_pdb_fragment(args.target_selection)
    model_selection = select_pdb_fragment(args.model_selection)

    if target_selection:
        if args.verbose:
            print('  target_selection: #', args.target_selection, target_selection)
            ts = target_selection
            print(ts)
            resides_in_total = 0
            for i in target_selection:
                print((i, len(ts[i]))) # chain string, a list of residues
                resides_in_total += len(ts[i])
            print('in total:', resides_in_total)

    if model_selection:
        if args.verbose:
            print('  model_selection:  ', len(model_selection), args.model_selection, model_selection)
            resides_in_total = 0
            for i in model_selection:
                print(i, len(model_selection[i])) # chain string, a list of residues
                resides_in_total += len(model_selection[i])
            print('in total:', resides_in_total)

    if args.target_ignore_selection:
        target_ignore_selection = select_pdb_fragment_pymol_style(args.target_ignore_selection, args.verbose)
    else:
        target_ignore_selection = None

    if args.model_ignore_selection:
        model_ignore_selection = select_pdb_fragment_pymol_style(args.model_ignore_selection, args.verbose)
    else:
        model_ignore_selection = None

    if  args.target_ignore_selection:
        if args.verbose: print('  target_ignore_selection: ', args.target_ignore_selection)
    if  args.model_ignore_selection:
        if args.verbose: print('  model_ignore_selection:  ', args.model_ignore_selection)

    models = get_rna_models_from_dir(input_files)

    print('# of models:', len(models))

    f = open(rmsds_fn, 'w')
    #t = 'target:' + os.path.basename(target_fn) + ' , rmsd_all\n'
    if args.target_column_name:
        t = 'fn,' + os.path.basename(args.target_fn) + '\n'
    else:
        t = 'fn,rmsd_all\n'

    c = 1
    for r1 in models:
        if method == 'align' or method == 'fit':
            rmsd_curr, atoms = calc_rmsd_pymol(r1, target_fn, method, args.verbose)
        else:
            rmsd_curr, atoms = calc_rmsd(r1, target_fn, target_selection, target_ignore_selection, model_selection, model_ignore_selection, args.verbose)
        r1_basename = os.path.basename(r1)
        if args.print_progress: print(r1_basename, rmsd_curr, atoms)
        t += r1_basename + ',' + str(round(rmsd_curr,3)) + ' '
        c += 1
        t += '\n'

    f.write(t)
    f.close()

    #print t.strip() # matrix

    print('# of atoms used:', atoms)

    df = pd.read_csv(rmsds_fn)
    df = df.round(2)
    if args.sort_results:
        df = df.sort_values('rmsd_all', ascending=True)
    if args.print_results:
        print(df)
    df.to_csv(rmsds_fn, sep=',', index=False)  # easy to set \t here!

    print('csv was created! ', rmsds_fn)
