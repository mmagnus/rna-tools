#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

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
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    if args.hires:
        headers = 'id,fn,' + ','.join(['farna_atom_pair_constraint_hires', 'farna_score_hires', 'farna_hbond_sr_bb_sc', 'farna_linear_chainbreak_hires', 'farna_fa_intra_rep', 'farna_rna_sugar_close', 'farna_hbond_sc', 'farna_lk_nonpolar', 'farna_fa_atr', 'farna_fa_elec_rna_phos_phos', 'farna_ch_bond', 'farna_fa_rep', 'farna_hbond_lr_bb_sc', 'farna_geom_sol', 'farna_rna_torsion'])
    else:
        headers = 'id,fn,' + ','.join(['farna_rna_base_axis', 'farna_rna_backbone_backbone', 'farna_rna_base_stack_axis', 'farna_rna_base_stagger', 'farna_rna_base_stack', 'farna_rna_base_pair', 'farna_rna_repulsive', 'farna_rna_vdw', 'farna_rna_base_backbone', 'farna_score_lowres', 'farna_rna_data_backbone', 'farna_linear_chainbreak', 'farna_rna_rg', 'farna_atom_pair_constraint'])
    print(headers)
    for i, f in enumerate(args.file):
        # mini false
        # print(f)
        farna = FARNA()
        try:
            result = farna.run(f, args.hires, args.verbose)  # False or True for min
        except:
            result = ''
        print(len(headers.split(',')) - 2 == len(result))
        # assert len(headers.split(',')) - 2 == len(result), 'error at ' + f
        print(','.join([str(i + 1), os.path.basename(f)]) + ',' + ','.join([str(x) for x in result]))
        #except:
        #   print(','.join([str(i + 1), os.path.basename(f)]) + ',' + 'error')
