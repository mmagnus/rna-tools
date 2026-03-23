#!/usr/bin/env python

"""ClaRNA app — pure-Python interface for INF calculation.

No subprocess calls. Everything runs in-process so it works anywhere
rna-tools is importable (Colab, OpenRNAFold, etc.).

Quick usage::

    from rna_tools.tools.clarna_app.rna_clarna_app import calc_inf

    scores = calc_inf('target.pdb', 'pred.pdb')
    print(scores['inf_all'])   # 0.706
    print(scores['inf_WC'])    # 0.865
"""

import argparse
import sys
import os
import tempfile

from rna_tools.tools.clarna_play.rna_clarna_run import run_clarna_direct
from rna_tools.tools.clarna_play.rna_clarna_compare import compare_clarna_direct
from rna_tools.tools.rna_convert_pseudoknot_formats.rna_pk_simrna_to_one_line import get_one_line


def clarna_run(fn, force=True, stacking=True, verbose=False):
    """Run ClaRNA annotation on a single PDB file (pure Python, no subprocess).

    Args:
        fn (str): path to PDB file
        force (bool): re-run even if .outCR already exists
        stacking (bool): include stacking interactions

    Returns:
        str: path to the generated .outCR file
    """
    if verbose:
        print('stacking', stacking)
    fn_out = fn + '.outCR'
    if os.path.isfile(fn_out) and not force:
        return fn_out

    clarna_opts = 'bp+stack' if stacking else 'bps'
    if verbose:
        print('clarna_run: direct call, opts=%s' % clarna_opts)
    result = run_clarna_direct(fn, clarna_opts=clarna_opts)
    with open(fn_out, 'w') as f:
        f.write(result)

    if os.stat(fn_out).st_size == 0:
        result = run_clarna_direct(fn, clarna_opts='bp+stack')
        with open(fn_out, 'w') as f:
            f.write(result)
    return fn_out


def clarna_compare(target_cl_fn, i_cl_fn, verbose=False):
    """Compare two .outCR files and return the raw score line (pure Python).

    Args:
        target_cl_fn (str): path to reference .outCR
        i_cl_fn (str): path to model .outCR

    Returns:
        str: space-separated line with filenames and 8 scores
    """
    scores = compare_clarna_direct(target_cl_fn, i_cl_fn)
    return scores['raw_line']


def calc_inf(target_pdb, model_pdb, stacking=True, force=True, verbose=False):
    """Calculate INF scores between two PDB structures.

    Pure-Python, no subprocess. Classifier libraries are cached after
    the first call, so repeated calls are fast.

    Args:
        target_pdb (str): path to reference PDB file
        model_pdb (str): path to predicted PDB file
        stacking (bool): include stacking interactions (default True)
        force (bool): re-run ClaRNA even if .outCR exists (default True)
        verbose (bool): print debug info

    Returns:
        dict: scores with keys inf_all, inf_stack, inf_WC, inf_nWC,
              sns_WC, ppv_WC, sns_nWC, ppv_nWC.
              Values are float or None (for NA).

    Example::

        from rna_tools.tools.clarna_app.rna_clarna_app import calc_inf
        scores = calc_inf('target.pdb', 'pred.pdb')
        print(scores['inf_all'])   # 0.706
        print(scores['inf_WC'])    # 0.865
    """
    target_cl_fn = clarna_run(target_pdb, force=force, stacking=stacking, verbose=verbose)
    model_cl_fn = clarna_run(model_pdb, force=force, stacking=stacking, verbose=verbose)
    return compare_clarna_direct(target_cl_fn, model_cl_fn)


def get_ClaRNA_output_from_dot_bracket(ss, temp=True, verbose=False):
    """Get ClaRNA .outCR output from a dot-bracket secondary structure string.

    Args:
        ss (str): secondary structure, optionally prefixed with chain like "A:((...))"

    Returns:
        str: path to generated .outCR file
    """
    from rna_tools.SecondaryStructure import parse_vienna_to_pairs

    if ss.find(':') > -1:
        chain, ss = ss.split(':')
    else:
        chain = 'A'

    pairs, pairs_pk = parse_vienna_to_pairs(ss, remove_gaps_in_ss=False)
    pairs += pairs_pk

    txt = 'Classifier: Clarna\n'
    txt += 'chains:  A 1 ' + str(len(ss)) + '\n'
    for bp in pairs:
        txt += '%s    %i   %s   %i          bp G C                  WW_cis   1 \n' % (chain, bp[0], chain, bp[1])
    if verbose:
        print(txt.strip())

    if temp:
        f = tempfile.NamedTemporaryFile()
        name = f.name
    else:
        name = 'target'

    foutCR = name + '.pdb.outCR'
    if verbose:
        print(foutCR)
    with open(foutCR, 'w') as ft:
        ft.write(txt)
    return foutCR


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', help="files", nargs='+')
    parser.add_argument('-f', "--force", dest="force", action="store_true",
                        help="force to run ClaRNA")
    parser.add_argument('-v', "--verbose", dest="verbose", action="store_true",
                        help="verbose")
    return parser


# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print(parser.print_help())
        sys.exit(1)

    for f in args.files:
        print(f)
        fn_out = clarna_run(f, args.force)
        # get_dot_bracket_from_ClaRNAoutput requires external ClaRNAwd_to_vienaSS binary
        import subprocess
        ClaRNA_play_path = os.path.dirname(os.path.abspath(__file__)) + '/../clarna_play'
        cmd = ClaRNA_play_path + '/lib/ClaRNAwd_to_vienaSS/ClaRNAwd_output_parser_get_SS ' + fn_out
        o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        mdb = o.stdout.read().strip()
        db = get_one_line(mdb.split('\n'))
        print(db)
