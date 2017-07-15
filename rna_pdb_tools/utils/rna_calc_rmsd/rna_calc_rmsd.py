#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
If you still have problem with various number of atoms, check out this issue: get_rnapuzzle_ready: how to deal with `Alternate location indicator (https://github.com/mmagnus/rna-pdb-tools/issues/30).

The program is using (https://github.com/charnley/rmsd).
"""
from __future__ import print_function

from rna_pdb_tools.utils.rna_calc_rmsd.lib.rmsd.calculate_rmsd import *
import sys
from rna_pdb_tools.utils.extra_functions.select_fragment import select_pdb_fragment_pymol_style, select_pdb_fragment
import optparse
import sys
import math
import glob
import re
import os

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
    
    if atomsQ != atomsP:
        sys.exit('Error: # of atoms is not equal target (' + b + '):' + str(atomsQ) + ' vs model (' + a + '):' + str(atomsP))
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

# main
if __name__ == '__main__':
    print('rmsd_calc_rmsd_to_target')
    print('-' * 80)
    
    optparser=optparse.OptionParser(usage="%prog [<options>] <pdb files (test_data/*)>")

    optparser.add_option('-t',"--target_fn", type="string",
                         dest="target_fn",
                         default='',
                         help="pdb file")

    optparser.add_option('',"--target_selection", type="string",
                         dest="target_selection",
                         default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")

    optparser.add_option('',"--target_ignore_selection", type="string",
                         dest="target_ignore_selection",
                         default='',
                         help="A/10/O2\'")
    
    optparser.add_option('',"--model_selection", type="string",
                         dest="model_selection",
                         default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")

    optparser.add_option('',"--model_ignore_selection", type="string",
                         dest="model_ignore_selection",
                         default='',
                         help="A/10/O2\'")
    
    optparser.add_option('-m',"--method", type="string",
                         dest="method",
                         default='all-atom-built-in',
                         help="align, fit")

    optparser.add_option('-o',"--rmsds_fn", type="string",
                         dest="rmsds_fn",
                         default='rmsds.csv',
                         help="ouput, matrix")
    
    optparser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="verbose")

    (opts, args)=optparser.parse_args()

    if len(sys.argv) == 1:
        print(optparser.format_help()) #prints help if no arguments
        sys.exit(1)

    input_files = args[:] # opts.input_dir
    rmsds_fn = opts.rmsds_fn
    target_fn = opts.target_fn
    method = opts.method
    print('method:', method)
    target_selection = select_pdb_fragment(opts.target_selection)
    model_selection = select_pdb_fragment(opts.model_selection)
    if target_selection:
        if opts.verbose:
            print('  target_selection: #', opts.target_selection, target_selection)
            ts = target_selection
            print(ts)
            resides_in_total = 0
            for i in target_selection:
                print((i, len(ts[i]))) # chain string, a list of residues
                resides_in_total += len(ts[i])
            print('in total:', resides_in_total)
    if model_selection:
        if opts.verbose:
            print('  model_selection:  ', len(model_selection), opts.model_selection, model_selection)
            resides_in_total = 0
            for i in model_selection:
                print(i, len(model_selection[i])) # chain string, a list of residues
                resides_in_total += len(model_selection[i])
            print('in total:', resides_in_total)
    if opts.target_ignore_selection:
        target_ignore_selection = select_pdb_fragment_pymol_style(opts.target_ignore_selection)
    else:
        target_ignore_selection = None
        
    if opts.model_ignore_selection:
        model_ignore_selection = select_pdb_fragment_pymol_style(opts.model_ignore_selection)
    else:
        model_ignore_selection = None

    if  opts.target_ignore_selection:
        if opts.verbose: print('  target_ignore_selection: ', opts.target_ignore_selection)
    if  opts.model_ignore_selection:
        if opts.verbose: print('  model_ignore_selection:  ', opts.model_ignore_selection)

    models = get_rna_models_from_dir(input_files)        

    print('# of models:', len(models))

    f = open(rmsds_fn, 'w')
    #t = 'target:' + os.path.basename(target_fn) + ' , rmsd_all\n'
    t = 'fn,rmsd_all\n'

    c = 1
    for r1 in models:
        if method == 'align' or method == 'fit':
            rmsd_curr, atoms = calc_rmsd_pymol(r1, target_fn, method)
        else:
            rmsd_curr, atoms = calc_rmsd(r1, target_fn, target_selection, target_ignore_selection, model_selection, model_ignore_selection, opts.verbose)
        r1_basename = os.path.basename(r1)
        print(r1_basename, rmsd_curr, atoms)
        t += r1_basename + ',' + str(round(rmsd_curr,3)) + ' '
        c += 1
        t += '\n'
            
    f.write(t)
    f.close()

    #print t.strip() # matrix

    print('# of atoms used:', atoms)
    if opts.rmsds_fn:
        print('csv was created! ', rmsds_fn)
