#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rmsd_calc_trafl - calculate RMSD of transition A->B based on a SimRNA trajectory

After this script, run::

    rna_cal_rmsd_trafl_plot.py rmsd.txt

to get a plot like this:

.. image:: ../pngs/rmsd_transition.png

Prepare structures::

     $ SimRNA -p 17_Das_2_rpr.pdb -n 0 -o 17_Das_2_rpr_n0 # no trafl, trafl will be added
     $ SimRNA -p 5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped.pdb -n 0 -o 5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped
     #(struc must be (~CG~) nope. It has to be a trajectory!)

and run::

     $ rmsd_calc_trafl.py 17_Das_2_rpr.pdb.trafl 17_Das_2_rpr_n0.trafl 5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped_n0.trafl rp17_rmsd.txt
     > calc_rmsd_to_1st_frame
       /Users/magnus/work/opt/simrna/SimRNA_64bitIntel_MacOSX_staticLibs/calc_rmsd_to_1st_frame 17_Das_2_rpr.pdb.trafl 17_Das_2_rpr.pdb_rmsd_e
     < rmsd_out: 17_Das_2_rpr.pdb_rmsd_e
     > struc: 17_Das_2_rpr_n0.trafl 2
     > trafl: 17_Das_2_rpr.pdb.trafl 48
     % saved: 17_Das_2_rpr.pdb.trafl_17_Das_2_rpr_n0.trafl
     > calc_rmsd_to_1st_frame
       /Users/magnus/work/opt/simrna/SimRNA_64bitIntel_MacOSX_staticLibs/calc_rmsd_to_1st_frame 17_Das_2_rpr.pdb.trafl_17_Das_2_rpr_n0.trafl 17_Das_2_rpr.pdb_rmsd_e_17_Das_2_rpr_n0_rmsd_e
     < rmsd_out: 17_Das_2_rpr.pdb_rmsd_e_17_Das_2_rpr_n0_rmsd_e
     > struc: 5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped_n0.trafl 2
     > trafl: 17_Das_2_rpr.pdb.trafl 48
     % saved: 17_Das_2_rpr.pdb.trafl_5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped_n0.trafl
     > calc_rmsd_to_1st_frame
       /Users/magnus/work/opt/simrna/SimRNA_64bitIntel_MacOSX_staticLibs/calc_rmsd_to_1st_frame 17_Das_2_rpr.pdb.trafl_5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped_n0.trafl 17_Das_2_rpr.pdb_rmsd_e_5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped_n0_rmsd_e
     < rmsd_out: 17_Das_2_rpr.pdb_rmsd_e_5k7c_clean_onechain_renumber_as_puzzle_rpr_rmGapped_n0_rmsd_e
      0.000 -695.634
      0.000 -551.093
    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     < out: rp17_rmsd.txt

.. warning:: calc_rmsd_to_1st_frame (SimRNA) is required and the path to the binary file is defined in ``config_local``.
"""
from __future__ import print_function
import argparse
import sys
import os
import itertools
import subprocess


class ExceptionRmsdCalcTrafl(Exception):
    pass


# Settings
try:
    from config_local import calc_rmsd_to_1st_frame_exec
except ImportError:
    print("""Create config_local file in the folder of this tool, and set the path to calc_rmsd_to_1st_frame_exec
e.g. calc_rmsd_to_1st_frame_exec="/Users/magnus/work/opt/simrna/SimRNA_64bitIntel_MacOSX_staticLibs/calc_rmsd_to_1st_frame""")


def add_to_trafl(trafl, struc_trafl):
    """Take trafl and struc_trafl and merge it into trafl (out_trafl).

    Args:

         trafl: trafl filename
         struc_trafl: struc_trafl filename (SimRNA trafl file of one structure)

    Returns:

         str: out_trafl, a path to an output file
    """
    print('\nadd_to_trafl')
    struc = open(struc_trafl).read().strip()
    # print struc
    # print struc
    nr_of_files_structure = len(struc.split('\n'))
    assert(nr_of_files_structure == 2,
           'Structure should be in SimRNA trajectory format (with only one frame)')
    print(' input struc: and # of frames (should be only 2!)', struc_trafl, nr_of_files_structure)
    trafl_txt = open(trafl).read().strip()
    print(' add the structure to this trafl: and # of frames', trafl, len(trafl_txt.split('\n')))

    # if len(struc.split('\n')) != len(trafl_txt.split('\n')):
    #    raise Exception('# of atoms in structure != # of atoms in trafl')

    # out_trafl = '/archive/magnus/rmsd_calc_trafl_tmp/' + os.path.basename(trafl) + "_" + os.path.basename(struc_trafl) # temp hack
    # out_trafl = '/archive/magnus/rmsd_calc_trafl_tmp/' + os.path.basename(trafl) + "_" + os.path.basename(struc_trafl) # temp hack
    out_trafl = trafl + "_" + os.path.basename(struc_trafl)
    f = open(out_trafl, 'w')
    f.write(struc + '\n' + trafl_txt)
    f.close()

    print(' saved:', out_trafl)
    return out_trafl


def calc_rmsd_to_1st_frame(trafl):
    """Calc rmsd to the 1st frame of a trafl.

    :param trafl: trafl filename, string
    """
    print('\ncalc_rmsd_to_1st_frame')
    rmsd_out = trafl.replace('.trafl', '_rmsd_e')
    cmd = calc_rmsd_to_1st_frame_exec + ' ' + trafl + ' ' + rmsd_out
    print(' calc_rmsd_to_1st_frame')
    print('  ', cmd)

    o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stdout.flush())  # subprocess.PIPE)
    out = o.stdout.read().strip()
    #err = o.stderr.read().strip()
    # if 'warning' in err:
    #    print('!!! error in reading a trafl', err, out)
    #    sys.exit(0)
    print(' rmsd_out:', rmsd_out)
    return rmsd_out


def head_trafl(trafl, n):
    """ n number of structures, so an output file will have n x 2 lines"""
    fi = open(trafl, 'r')
    trafl_out = trafl.replace('.trafl', '.head' + str(n) + '.trafl')
    fo = open(trafl_out, 'w')
    n = n * 2
    c = 0
    for l in fi:
        fo.write(l)
        c += 1
        if c == n:
            break
    fo.close()
    print(' < saved:', trafl_out)
    return trafl_out


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('trafl', help="trafil")
    parser.add_argument('struc1', help="structure A")
    parser.add_argument('struc2', help="structure B")
    parser.add_argument('rmsds_fn', help="output file")
    return parser


# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    # if you want to head, then uncomment
    #trafl = head_trafl(trafl, 10000000000000000000000)
    calc_rmsd_to_1st_frame(args.trafl)

    trafl = args.trafl
    struc1 = args.struc1
    struc2 = args.struc2
    rmsds_fn = args.rmsds_fn

    # struc 1
    out_trafl = add_to_trafl(trafl, struc1)
    rmsd1_fn = calc_rmsd_to_1st_frame(out_trafl)

    # struc 2
    out_trafl2 = add_to_trafl(trafl, struc2)
    rmsd2_fn = calc_rmsd_to_1st_frame(out_trafl2)

    rmsd1 = open(rmsd1_fn).read().split('\n')
    rmsd2 = open(rmsd2_fn).read().split('\n')

    print(rmsd1[0])
    print(rmsd2[0])

    print('x' * 80)

    c = 0

    fo = open(rmsds_fn, 'w')

    # header
    # fo.write('smk_ready_samsoaked_chainB' + '\t' + 'smk_ready_samsoaked_chainA+SAM' + '\n') #th.basename(struc1).replace(' + '\t' + os.path.basename(struc2) + '\n')
    fo.write(os.path.basename(struc1) + '\t' + os.path.basename(struc2) + '\n')

    for item in itertools.izip(rmsd1, rmsd2):
        x = str(item).split()  # ('  0.000 -502.973', '  0.000 -496.943')
        # print x #["('", '0.000', "-502.973',", "'", '0.000', "-496.943')"]
        if x[0] == "('',":
            break
        r1 = float(x[1])
        r2 = float(x[4])
        # print r1, r2
        if r1 == 0 and r2 == 0:  # skip 0 0
            continue
        fo.write(str(r1) + '\t' + str(r2) + '\n')
        c += 1
    fo.close()
    print(' < out:', rmsds_fn)
