#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example::

    rna_mq_dfire.py AA/*.pdb
    id,fn,dfire
    1,1msy_A-000017_AA.pdb,-16181.328803
    2,1msy_A-000176_AA.pdb,-18172.496594
    3,1msy_A-000205_AA.pdb,-15526.172071

Output will be printed to stdout.
"""
from __future__ import print_function
import argparse
from rna_tools.tools.mq.SimRNA.SimRNA import SimRNA
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument('-', "--", help="", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    print('id,fn,' + ','.join(['simrna_steps', 'simrna_total_energy', 'simrna_base_base', 'simrna_short_stacking', 'simrna_base_backbone',  'simrna_local_geometry', 'simrna_bonds_dist_cp', 'simrna_bonds_dist_pc', 'simrna_flat_angles_cpc', 'simrna_flat_angles_pcp', 'simrna_tors_eta_theta', 'simrna_sphere_penalty', 'simrna_chain_energy']))
    for i, f in enumerate( args.file):
        Sim = SimRNA()
        #result = Sim.run('test' + os.sep + 'decoy3308.pdb', '1')#32000')
        result = Sim.run(f, '0', verbose=args.verbose) # two chains #32000')
        print(','.join([str(i + 1), os.path.basename(f)]) + ',' + ','.join(result))
