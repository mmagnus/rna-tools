#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from rna_tools.tools.mq.FARNA.FARNA import FARNA
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
        print('id,fn,farna_score_hires,farna_fa_atr,farna_fa_rep,farna_fa_intra_rep,farna_lk_nonpolar,farna_fa_elec_rna_phos_phos,farna_rna_torsion,farna_suiteness_bonus,farna_rna_sugar_close,farna_fa_stack,farna_stack_elec,farna_geom_sol_fast,farna_bond_sr_bb_sc,farna_hbond_lr_bb_sc,farna_hbond_sc,farna_ref,farna_free_suite,farna_free_2HOprime,farna_intermol,farna_other_pose,farna_loop_close,farna_linear_chainbreak_hires')
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
        farna = FARNA()
        result = farna.run(f, args.hires, args.verbose)  # False or True for min
        print(','.join([str(i + 1), os.path.basename(f)]) + ',' + ','.join(result))
