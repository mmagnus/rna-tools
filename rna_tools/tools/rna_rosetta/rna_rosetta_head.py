#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_roseta_head.py - get a head of a Rosetta silent file.

Example::

     $ rna_rosetta_head.py -n 10000 thf.out
     # a new file will be created, thf_10000.out

Silent file::

    [peyote2] thf head -n 100 thf.out
    SEQUENCE: ggagaguagaugauucgcguuaagugugugugaaugggaugucgucacacaacgaagcgagagcgcggugaaucauugcauccgcucca
    SCORE:     score    rna_data_backbone    rna_vdw    rna_base_backbone    rna_backbone_backbone    rna_repulsive    rna_base_pair    rna_base_axis    rna_base_stagger    rna_base_stack    rna_base_stack_axis     rna_rg    atom_pair_constraint    linear_chainbreak       N_WC      N_NWC       N_BS description
    REMARK BINARY SILENTFILE
    SCORE:  -601.975                0.000     31.113              -16.960                   -3.888           20.501         -302.742          -38.531            -158.004           -80.764               -110.053     23.750                   0.000               33.601         32          6         86    S_000001_000
    FOLD_TREE  EDGE 1 4 -1  JEDGE 4 85 1  C4   C2   END  EDGE 4 5 -1  EDGE 85 80 -1  EDGE 85 89 -1  JEDGE 80 40 5  C4   C2   END  EDGE 80 78 -1  EDGE 40 43 -1  EDGE 40 33 -1  JEDGE 33 45 4  C4   C2   END  EDGE 45 54 -1  EDGE 45 44 -1  EDGE 33 20 -1  JEDGE 20 65 3  C2   C4   END  EDGE 65 67 -1  EDGE 20 17 -1  EDGE 65 63 -1  JEDGE 17 69 2  C4   C2   END  JEDGE 63 58 6  C4   C2   END  EDGE 17 6 -1  EDGE 58 59 -1  EDGE 63 60 -1  EDGE 69 77 -1  EDGE 58 55 -1  EDGE 69 68 -1  S_000001_000
    RT -0.987743 0.139354 0.0703103 0.135963 0.989404 -0.0509304 -0.0766626 -0.0407466 -0.996224 6.25631 0.103544 0.0647696   S_000001_000
    RT -0.98312 0.1587 -0.091045 0.166923 0.981743 -0.0912024 0.074909 -0.10486 -0.991662 5.89962 -1.95819 -0.1075   S_000001_000
    RT -0.987645 0.154078 0.0285994 0.153854 0.988044 -0.00989514 -0.0297821 -0.00537275 -0.999542 6.13138 1.047 0.115722   S_000001_000
    RT -0.989884 0.140722 0.0180554 0.140036 0.989532 -0.0348618 -0.0227723 -0.0319807 -0.999229 6.21787 0.201588 -0.0264223   S_000001_000
    RT -0.988455 0.134013 0.0706914 0.134924 0.990822 0.00825457 -0.0689364 0.0176973 -0.997464 6.19447 0.189237 0.125791   S_000001_000
    RT -0.990412 0.138036 0.00546261 0.137927 0.990299 -0.0168644 -0.00773751 -0.0159492 -0.999843 6.25456 0.0842849 -0.0135807   S_000001_000
    ANNOTATED_SEQUENCE: g[RGU:LowerRNA:Virtual_Phosphate]gaga[RAD:rna_cutpoint_lower]g[RGU:rna_cutpoint_upper]uagaugauucgcguuaagugugugugaaugggauguc[RCY:rna_cutpoint_lower]g[RGU:rna_cutpoint_upper]ucacacaacg[RGU:rna_cutpoint_lower]a[RAD:rna_cutpoint_upper]agcg[RGU:rna_cutpoint_lower]a[RAD:rna_cutpoint_upper]gagcgcg[RGU:rna_cutpoint_lower]g[RGU:rna_cutpoint_upper]ugaaucauu[URA:rna_cutpoint_lower]g[RGU:rna_cutpoint_upper]cauccgcucca[RAD:UpperRNA] S_000001_000
    Le3nY9smsa4zEMGdvAA+z+e

It seems to work::

    -rw-rw-r--   1 magnus users 474M 2017-08-06 05:25 thf_10000.out
    -rw -rw-r--   1 magnus users 947M 2017-08-06 04:54 thf.out

    [peyote2] thf rna_rosetta_n.py thf_10000.out
    10000

"""
from __future__ import print_function
import subprocess
import argparse
import sys


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-n', '--nstruc', type=int, default=10000)
    parser.add_argument('file', help='ade.out')
    return parser


def run():
    """Pipline for modeling RNA."""
    args = get_parser().parse_args()
    nfile = ''
    c = -2
    for l in open(args.file):
        # print(l)
        if l.startswith('SCORE'):
            c += 1
            if c == int(args.nstruc):
                f = open(args.file.replace('.out', '_' + str(args.nstruc) + '.out'), 'w')
                # print(nfile)
                f.write(nfile)
                f.close()
                sys.exit(1)
        nfile += l


# main
if __name__ == '__main__':
    run()
