#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

    rna_mq_farfar2.py -q *.pdb | tee ff2.csv


    (base) ➜  NNN parallel rna_mq_farfar2.py -q  {} ::: *.pdb | tee ff2.csv
    1,nnn.out.6.pdb,-176.023,10.467,-9.662,-6.419,18.463,-122.033,-13.794,-40.49,-12.481,-15.034,14.96,0.0,0.0
    1,nnn.out.1.pdb,-160.603,11.138,-9.657,-6.413,17.474,-113.153,-12.704,-35.709,-11.96,-14.589,14.972,0.0,0.0
    1,nnn.out.3.pdb,-169.382,11.36,-10.322,-6.172,17.853,-119.688,-13.435,-39.454,-10.96,-13.533,14.969,0.0,0.0
    1,nnn.out.2.pdb,-164.123,10.651,-9.19,-6.41,17.707,-116.158,-13.07,-35.984,-11.96,-14.676,14.965,0.0,0.0
    1,nnn.out.4.pdb,-162.867,11.561,-10.142,-6.459,17.957,-116.398,-13.05,-35.971,-11.521,-13.808,14.963,0.0,0.0
    1,nnn.out.10.pdb,-172.812,11.311,-10.829,-4.577,18.429,-121.431,-13.678,-38.904,-12.617,-15.459,14.943,0.0,0.0
    1,nnn.out.5.pdb,-168.892,11.228,-10.404,-6.034,17.317,-118.941,-13.372,-38.279,-11.498,-13.875,14.967,0.0,0.0
    1,nnn.out.7.pdb,-154.926,12.819,-10.282,-6.45,17.539,-110.714,-12.364,-35.907,-10.96,-13.579,14.972,0.0,0.0
    1,nnn.out.9.pdb,-170.702,11.72,-10.09,-6.449,17.451,-122.965,-13.779,-37.119,-10.96,-13.479,14.967,0.0,0.0
    1,nnn.out.8.pdb,-162.495,11.094,-10.431,-5.457,18.034,-115.686,-12.956,-36.254,-11.641,-14.151,14.952,0.0,0.0
    1,pre2nd-delt2-3.pdb,-156.159,22.102,-7.019,-4.175,16.939,-117.089,-12.63,-41.987,-11.993,-15.285,14.978,0.0,0.0

Header for hires:

id,fn,ff2_score_hires,ff2_fa_atr,ff2_fa_rep,ff2_fa_intra_rep,ff2_lk_nonpolar,ff2_fa_elec_rna_phos_phos,ff2_rna_torsion,ff2_suiteness_bonus,ff2_rna_sugar_close,ff2_fa_stack,ff2_stack_elec,ff2_geom_sol_fast,ff2_bond_sr_bb_sc,ff2_hbond_lr_bb_sc,ff2_hbond_sc,ff2_ref,ff2_free_suite,ff2_free_2HOprime,ff2_intermol,ff2_other_pose,ff2_loop_close,ff2_linear_chainbreak_hires

for lowres:

id,fn,farna_score_lowres,farna_rna_vdw,farna_rna_base_backbone,farna_rna_backbone_backbone,farna_rna_repulsive,farna_rna_base_pair,farna_rna_base_axis,farna_rna_base_stagger,farna_rna_base_stack,farna_rna_base_stack_axis,farna_rna_rg,farna_atom_pair_constraint,farna_linear_chainbreak

"""
from __future__ import print_function
import argparse
from rna_tools.tools.mq.FARFAR2.FARFAR2 import FARFAR2
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-q", "--quiet",
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

    if not args.quiet:
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
