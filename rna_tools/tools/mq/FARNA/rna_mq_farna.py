#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Scores for `test_data/1y26X_output4_01-000001_AA.pdb`::

    ------------------------------------------------------------
     Scores                       Weight   Raw Score Wghtd.Score
    ------------------------------------------------------------
     fa_atr                       0.230    -802.076    -184.477
     fa_rep                       0.120     105.478      12.657
     fa_intra_rep                 0.003     531.095       1.540
     lk_nonpolar                  0.320      23.335       7.467
     fa_elec_rna_phos_phos        1.050       1.643       1.725
     ch_bond                      0.420    -263.589    -110.707
     rna_torsion                  0.100      46.967       4.697
     rna_sugar_close              0.700      11.901       8.330
     hbond_sr_bb_sc               0.620      -1.012      -0.627
     hbond_lr_bb_sc               3.400      -2.254      -7.664
     hbond_sc                     3.400    -118.246    -402.037
     geom_sol                     0.620     332.855     206.370
     atom_pair_constraint         1.000       0.000       0.000
     linear_chainbreak            5.000       0.000       0.000
    ---------------------------------------------------
     Total weighted score:                     -462.726

    /Users/magnus/work/opt/rosetta/rosetta_bin_mac_2016.13.58602_bundle/rna_minimize.static.macosclangrelease -database /Users/magnus/work/opt/rosetta/rosetta_bin_mac_2016.13.58602_bundle/database    -ignore_zero_occupancy false  -s /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpaf0ic8g7/query.pdb -out:file:silent /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpaf0ic8g7/SCORE.out
    100% (1 of 1) |#########################################################################################################################################################| Elapsed Time: 0:03:37 ETA:  00:00:001y26X_output4_01-000001_AA.pdb {'FARNA_hires': ['-462.726', '-184.477', '12.657', '1.54', '7.467', '1.725', '-110.707', '4.697', '8.33', '-0.627', '-7.664', '-402.037', '206.37', '0.0', '0.0'], 'length': 

"""
from __future__ import print_function
import argparse
from rna_tools.tools.mq.FARNA.FARNA import FARNA
import os


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-r", "--hires",
                        action="store_true", help="")
    parser.add_argument("-s", "--system",
                        action="store_true", help="")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    if args.hires:
        headers = 'id,fn,' + ','.join(['farna_score_hires', 'farna_fa_atr', 'farna_fa_rep', 'farna_fa_intra_rep',
                                       'farna_lk_nonpolar',
                                       'farna_fa_elec_rna_phos_phos',
                                       'farna_ch_bond',
                                       'farna_rna_torsion',
                                       'farna_rna_sugar_close',
                                       'farna_hbond_sr_bb_sc',
                                       'farna_hbond_lr_bb_sc',
                                       'farna_hbond_sc',
                                       'farna_geom_sol',
                                       'farna_atom_pair_constraint_hires',
                                       'farna_linear_chainbreak_hires'])
    else:
        headers = 'id,fn,' + ','.join(['farna_score_lowres',
                                       'farna_rna_data_backbone',
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
                                       'farna_linear_chainbreak',
                                       ])
    print(headers)
    for i, f in enumerate(args.file):
        # mini false
        # print(f)
        farna = FARNA()
        try:
            result = farna.run(f, args.hires, args.verbose, args.system)  # False or True for min
        except:
            result = ''
        print(','.join([str(i + 1), os.path.basename(f)]) + ',' + ','.join([str(x) for x in result]))
