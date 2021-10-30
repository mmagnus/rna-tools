#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from rna_tools.tools.mq.FARFAR2.FARFAR2 import FARFAR2
import os
import pandas as pd

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('-d', "--done", help="a csv with already done scores", default="")
    parser.add_argument("-r", "--hires",
                        action="store_true", help="")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    if args.hires:
        print('id,fn,ff2_score_hires,ff2_fa_atr,ff2_fa_rep,ff2_fa_intra_rep,ff2_lk_nonpolar,ff2_fa_elec_rna_phos_phos,ff2_rna_torsion,ff2_suiteness_bonus,ff2_rna_sugar_close,ff2_fa_stack,ff2_stack_elec,ff2_geom_sol_fast,ff2_bond_sr_bb_sc,ff2_hbond_lr_bb_sc,ff2_hbond_sc,ff2_ref,ff2_free_suite,ff2_free_2HOprime,ff2_intermol,ff2_other_pose,ff2_loop_close,ff2_linear_chainbreak_hires')
    else:
        print('id,fn,' + ','.join(['farna_score_lowres',
             'farna_rna_vdw',
             'farna_rna_base_backbone',
             'farna_rna_backbone_backbone',
             'farna_rna_repulsive',
             'farna_rna_base_pair',
             'farna_rna_base_axis',
             'farna_rna_base_stagger',
             'farna_rna_base_stack',
             'farna_rna_base_stack_axis',
             'farna_rna_rg',
             'farna_atom_pair_constraint',
             'farna_linear_chainbreak']))
#                                   'farna_score', 'farna_rna_rg', 'farna_rna_vdw', 'farna_rna_base_backbone',               'farna_rna_backbone_backbone', 'farna_rna_repulsive', 'farna_rna_base_pair_pairwise',               'farna_rna_base_pair', 'farna_rna_base_axis', 'farna_rna_base_stagger',               'farna_rna_base_stack', 'farna_rna_base_stack_axis', 'farna_atom_pair_constraint']))
    for i, f in enumerate(args.file):
        # mini false
        # print(f)
        farna = FARFAR2()
        result = farna.run(f, args.hires, args.verbose)  # False or True for min
        print(','.join([str(i + 1), os.path.basename(f)]) + ',' + ','.join([str(x) for x in result]))
