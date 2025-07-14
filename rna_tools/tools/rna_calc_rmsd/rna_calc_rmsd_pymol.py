#-*- coding: utf-8 -*-
"""rna_calc_rmsd_pymol.py - calculate RMSDs of structures with PyMOL
"""
from __future__ import print_function

import warnings
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)

import sys
import argparse
import sys
import math
import glob
import re
import os

try:
    import __main__
    __main__.pymol_argv = ['pymol', '-qc']
    import pymol  # import cmd, finish_launching

    from pymol import cmd
    #cmd.feedback("disable")
    cmd.set("logging", 0)  # optional; turns off logging
    cmd.set("suspend_updates", 1)  # optional; disables GUI redraws (if in GUI)
    pymol.finish_launching()

except ImportError:
    print('calc_rmsd_pymol: you need to have installed PyMOL')
    sys.exit(0) # no error



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

def calc_rmsd_pymol(pdb1, pdb2, method, verbose=False):
    """Calculate rmsd using PyMOL. Two methods are available: align and fit

    See:

    -  Align: <http://www.pymolwiki.org/index.php/Align>
    -  Fit:   <http://www.pymolwiki.org/index.php/Fit>

    Align can return a list with 7 items:

    - RMSD after refinement
    - Number of aligned atoms after refinement
    - Number of refinement cycles
    - RMSD before refinement
    - Number of aligned atoms before refinement
    - Raw alignment score
    - Number of residues aligned

    in this version of function, the function returns `RMSD before refinement`.

    Install on OSX: ``brew install brewsci/bio/pymol`` or get 

    If you have a problem::

      Match-Error: unable to open matrix file '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/data/pymol/matrices/BLOSUM62'.

    then find BLOSUM62, e.g.::

        mdfind -name BLOSUM62 | grep pymol
        /Users/magnus/miniconda2/envs/py37/lib/python3.7/site-packages/pymol/pymol_path/data/pymol/matrices/BLOSUM62
        /usr/local/Cellar/pymol/2.4.0_3/libexec/lib/python3.9/site-packages/pymol/pymol_path/data/pymol/matrices/BLOSUM62
        /Users/magnus/miniconda2/pkgs/pymol-2.4.2-py37h06d7bae_0/share/pymol/data/pymol/matrices/BLOSUM62
        /Users/magnus/work/opt/pymol-open-source/data/pymol/matrices/BLOSUM62

    and then define ``PYMOL_DATA`` in your .bashrc/.zshrc, e.g.::

       export PYMOL_DATA="/Users/magnus/work/opt/pymol-open-source/data/pymol"

     """

    pymol.cmd.reinitialize()
    pymol.cmd.delete('all')
    pymol.cmd.load(pdb1, 's1')
    pymol.cmd.load(pdb2, 's2')

    if method == 'align':
        # experiments with align <https://pymolwiki.org/index.php/Align>
        # quiet = 0/1: suppress output {default: 0 in command mode, 1 in API}
        # (4.130036354064941, 60, 3, 4.813207626342773, 64, 30.0, 3)
        values = pymol.cmd.align('s1', 's2',quiet=1, object='aln')

        from Bio import pairwise2

        def get_sequence(obj):
            """Extract 1-letter FASTA sequence from PyMOL object"""
            fasta = cmd.get_fastastr(obj)
            lines = fasta.splitlines()
            return ''.join(lines[1:])  # Skip the >header line

        def calc_sequence_identity(seq1, seq2):
            """Calculate global sequence identity"""
            alignments = pairwise2.align.globalxx(seq1, seq2)
            best = alignments[0]
            matches = sum(a == b for a, b in zip(best.seqA, best.seqB))
            return matches / max(len(seq1), len(seq2))

        seq1 = get_sequence("s1")
        seq2 = get_sequence("s2")
        print(seq1, seq2)
        
        identity = calc_sequence_identity(seq1, seq2)
        print(f"Sequence identity: {identity:.2%}")

        print(values)
        return values[0], values[3], identity # (,#0) #, pymol.cmd.align('s1','s2')[4])
        #raw_aln = pymol.cmd.get_raw_alignment('aln')
        #print raw_aln
        #for idx1, idx2 in raw_aln:
        #    print '%s`%d -> %s`%d' % tuple(idx1 + idx2)
        #pymol.cmd.save('aln.aln', 'aln')

    if method == 'fit':
        return (pymol.cmd.fit('s1', 's2'), 'fit')

def calc_rmsd(a, b, target_selection, target_ignore_selection, model_selection, model_ignore_selection, way, verbose):
    """
    Calculate RMSD between two XYZ files

    by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Andersen Bratholm <larsbratholm@gmail.com>
    project: https://github.com/charnley/rmsd
    license: https://github.com/charnley/rmsd/blob/master/LICENSE

    a is model
    b is target

    :params: a = filename of structure a
    :params: b = filename of structure b

    :return: rmsd, number of atoms
    """
    if verbose: print('in:', a)

    atomsP, P, atoms = get_coordinates(a, model_selection, model_ignore_selection, 'pdb', True, way)
    atomsQ, Q, atoms = get_coordinates(b, target_selection,target_ignore_selection,  'pdb', True, way)

    if verbose:
        print(atomsP, P)
        print(atomsQ, Q)

    if atomsQ != atomsP:
        print('Error: number of atoms is not equal target (' + b + '):' + str(atomsQ) + ' vs model (' + a + '):' + str(atomsP))
        return (-1,0) # skip this RNA
    # Calculate 'dumb' RMSD
    normal_rmsd = rmsd(P, Q, atoms)
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

    return round(kabsch_rmsd(P, Q, atoms),2), atomsP

def get_parser():
    import argparse
    class SmartFormatter(argparse.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()  
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=SmartFormatter)#formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t',"--target-fn",
                         default='', required = True,
                         help="pdb file")

    parser.add_argument('--ignore-files', help='files to be ingored, .e.g, \'solution\'', default='')

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

    parser.add_argument('--way', help="""R|c1p = C1'
backbone = P OP1 OP2 O5' C5' C4' C3' O3'
po = P OP1 OP2
no-backbone = all - po
bases, backbone+sugar, sugar
pooo = P OP1 OP2 O5'
alpha = P OP1 OP2 O5' C5'
""", default='all')

    parser.add_argument("--name-rmsd-column", help="default: fn,rmsd, with this cols will be fn,<name-rmsd-column>")

    parser.add_argument("--target-column-name", action="store_true",
                        help="")

    parser.add_argument('files', help='files', nargs='+')

    return parser
# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    input_files = args.files  # opts.input_dir
    tmp = []
    if args.ignore_files:
        for f in input_files:
            if args.ignore_files in f:
                continue
            tmp.append(f)
        input_files = tmp

    rmsds_fn = args.rmsds_fn
    target_fn = args.target_fn
    target_fn_only = os.path.basename(target_fn)
    method = args.method
        
    if args.verbose:
        print('method:', method)

    models = get_rna_models_from_dir(input_files)

    if args.verbose:
        print('target:', target_fn)
        print('of models:', len(models))

    fl = open(rmsds_fn, 'w')
    #t = 'target:' + os.path.basename(target_fn) + ' , rmsd_all\n'
    if args.name_rmsd_column:
        t = 'fn,' + args.name_rmsd_column + '\n'
    elif args.target_column_name:
        t = 'fn,' + os.path.basename(args.target_fn) + '\n'
    else:
        t = 'fn,rmsd_all,seq_identity\n'

    if not args.verbose:
        t = ''
    c = 1
    for r1 in models:
        if method == 'align' or method == 'fit':
            rmsd_curr, atoms, seq_identity = calc_rmsd_pymol(r1, target_fn, method, args.verbose)
        r1_basename = os.path.basename(r1)
        if args.print_progress: print(r1_basename, rmsd_curr, atoms)
        line = f"{target_fn_only},{r1_basename},{str(round(rmsd_curr,3))},{str(seq_identity)}"
        #print(line)
        t += line
        c += 1
        t += '\n'

    fl.write(t)
    fl.close()

    if args.verbose:
        print('number of atoms used:', atoms)

    try:
        import pandas as pdx
        pd.set_option('display.max_rows', 1000)

    except:
        print(t.strip()) # matrix
        sys.exit(0)
        
    df = pd.read_csv(rmsds_fn)
    df = df.round(2)
    if args.sort_results:
        df = df.sort_values('rmsd_all', ascending=True)
    if args.print_results:
        print(df)
    df.to_csv(rmsds_fn, sep=',', index=False)  # easy to set \t here!

    # print('# csv was created! ', rmsds_fn)
