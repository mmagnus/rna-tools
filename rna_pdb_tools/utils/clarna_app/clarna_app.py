#!/usr/bin/env python

"""
usage::

  ./clarna_app.py ../../input/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb
  ../../input/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb
  ((((([[[[[[)))))........(.((....(]]]]]].)..(((......)))...)).)

Example:

 .. code-block:: python

    from rna_pdb_tools.utils.clarna_app import clarna_app
    if __name__ == '__main__':
        ss = '((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))'
        fnCRref = clarna_app.get_ClaRNA_output_from_dot_bracket(ss)
        f = '../rna_calc_rmsd/test_data/pistol/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb'
        fnCR = clarna_app.clarna_run(f, force=False)
        results = clarna_app.clarna_compare(fnCRref, fnCR)
        print results #
        #tmp_Z42i_..pdb.outCR     5k7c_clean_onechain_renumber_as_puzzle_srr.pdb.outCR      0.706      NA         0.865      NA         0.842      0.889      NA         0.000

.. warning:: Setup a bash variable: ClaRNA_play_path, and add ClaRNA_play to your $PATH (install ClaRNA_play https://gitlab.genesilico.pl/RNA/ClaRNA_play (internal GS gitlab server)"""

import argparse
import subprocess
import sys
import os
import tempfile
from rna_pdb_tools.utils.rna_convert_pseudoknot_formats.rna_pk_simrna_to_one_line import get_one_line

ClaRNA_play_path = os.environ['ClaRNA_play_path']


def clarna_run(fn, force):
    """Run ClaRNA run

    fn: str
     filename to analyze

    :return: a filename to ClaRNA output (fn + '.outCR')"""
    fn_out = fn + '.outCR'
    if os.path.isfile(fn_out) and not force:
        pass
    else:
        cmd = 'clarna_run.py -ipdb ' + fn + ' > ' + fn_out
        os.system(cmd)
    return fn_out

def get_dot_bracket_from_ClaRNAoutput(inCR, verbose):
    """In inCR file"""
    cmd = ClaRNA_play_path + '/lib/ClaRNAwd_to_vienaSS/ClaRNAwd_output_parser_get_SS ' + inCR
    o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std = o.stdout.read().strip()
    if verbose: std
    return std
    
def clarna_compare(target_cl_fn,i_cl_fn):
    """Run ClaRNA compare.

    :return: a list target, fn, scores

    Scores::

        inf_all      0.706
        inf_stack -999.999 -> NA
        inf_WC       0.865 
        inf_nWC   -999.999 -> NA
        SNS_WC       0.842
        PPV_WC       0.889
        SNS_nWC     NA
        PPV_nWC      0.000

    Example of the list::

       5k7c_clean_onechain_renumber_as_puzzle_srr.pdb     pistol_thrs0.50A_clust01-000001_AA.pdb
        0.642      NA         0.874      0.000      0.944      0.810      0.000      0.000s

    use ``results.split()[4]`` to get inf_WC"""

    cmd = 'clarna_compare.py -iref ' + target_cl_fn + ' -ichk ' + i_cl_fn
    o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std = o.stdout.read().strip()
    return std # std.replace('.outCR', '') 
    
def get_ClaRNA_output_from_dot_bracket(ss, temp=True, verbose=False):
    """
    Get dummy ClaRNA output out of dat bracket secondary structure (ss)

    :return: a filename to ClaRNA output"""
    from rna_pdb_tools.utils.secstruc.secstruc import ViennaStructure

    if ss.find(':') > -1:
        chain,ss = ss.split(':')
    else:
        chain = 'A'
        ss = ss

    v = ViennaStructure(ss)
    
    txt = 'Classifier: Clarna\n'
    txt += 'chains:  A 1 ' + str(len(ss)) + '\n'
    for bp in v.toPairs():
        txt += '%s    %i   %s   %i          bp G C                  WW_cis   1 \n' % (chain, bp[0], chain, bp[1])
    if verbose: print txt.strip()

    if temp:
        f = tempfile.NamedTemporaryFile()
        name = f.name
    else:
        name = '/tmp/target'

    foutCR = name + '.pdb.outCR'
    print  foutCR
    ft = open(foutCR, 'w')
    ft.write(txt)
    ft.close()
    return foutCR
    
def get_parser():
    parser =  argparse.ArgumentParser()#usage="%prog [<options>] <pdb files (test_data/*)>")
    parser.add_argument('files', help="files", nargs='+')
    parser.add_argument('-f',"--force",
                         dest="force",
                         action="store_true",
                        help="force to run ClaRNA")
    parser.add_argument('-v',"--verbose",
                         dest="verbose",
                         action="store_true",
                        help="verbose")

    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print parser.print_help()
        sys.exit(1)
        
    for f in args.files:
        print f
        fn_out = clarna_run(f, args.force)
        mdb = get_dot_bracket_from_ClaRNAoutput(fn_out, args.verbose) # multie line dot bracket
        db = get_one_line(mdb.split('\n'))
        print db
