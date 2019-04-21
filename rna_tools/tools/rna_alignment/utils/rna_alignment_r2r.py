#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_alignment_r2r.py - a simple wrapper for r2r
R2R
Software to speed the depiction of aesthetic consensus RNA secondary structures
Zasha Weinberg

Z. Weinberg and R. R. Breaker. R2R—software to speed the depiction of aesthetic consensus RNA secondary structures. BMC Bioinformatics, 21:3, 2011.

Example::

    rna_alignment_r2r.py Intron_gpII-A_BactA.stk

What's next: the script will be integrated into RNAalignment class
"""

from __future__ import print_function
import rna_tools.tools.rna_alignment.rna_alignment as ra
import argparse
import subprocess


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="fasta seq")
    parser.add_argument("--svg", default=True,
                        action="store_true", help="output also a svg file")
    parser.add_argument("-v", "--verbose", default=True,
                        action="store_true", help="be verbose")
    return parser


def exe(cmd, verbose):
    if verbose:
        print(cmd)
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


def r2r(alignment, svg, verbose):
    cons =  alignment + '.cons'
    pdf =  alignment.replace('.stk', '').replace('.sto', '') + '.pdf'
    svg =  alignment.replace('.stk', '').replace('.sto', '') + '.svg'
    cmd = 'r2r --GSC-weighted-consensus ' + alignment + ' ' + cons + ' 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1'
    out, err = exe(cmd, verbose)
    if verbose: print(out)
    if err:
        print(err)

    cmd = 'r2r ' + cons + ' ' + pdf
    out, err = exe(cmd, verbose)
    if err:
        print(err)
    if verbose: print(out)

    if svg:
        cmd = 'r2r ' + cons + ' ' + svg
        out, err = exe(cmd, verbose)
        if err:
            print(err)
        if verbose: print(out)


if __name__ == '__main__':
    # for debugging #
    #fn = 'test_data/RF00167.stockholm.sto'
    #ids = ['AL591975.1/251136-251218','CP000721.1/2204691-2204778']#ACCL02000010.1/116901-116991']
    parser = get_parser()
    args = parser.parse_args()

    r2r(args.file, args.svg, args.verbose)
