#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_simrna_extra.py - extract full atom structures from a SimRNA trajectory file

Options:

  SIMRNA_DATA_PATH has to be properly defined in ``rpt_config_local``.

"""
from rna_tools.tools.simrna_trajectory.simrna_trajectory import SimRNATrajectory
from rna_tools.rna_tools_config import SIMRNA_DATA_PATH

import argparse
import os

import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('rna_simrna_extract.log')
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)


def get_parser():
    """Get parser of arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-t', '--template',
        help="template PDB file used for reconstruction to full atom models", required=True)
    parser.add_argument('-f', '--trafl', help="SimRNA trafl file", required=True)
    parser.add_argument('-c', '--cleanup', action='store_true',
                        help="Keep only *_AA.pdb files, move *.ss_detected and *.pdb"
                        "to _<traj name folder>")
    parser.add_argument('-n', '--number_of_structures', help="", default=100)
    return parser


def get_data():
    """Get a link to SimRNA data folder in cwd."""
    cmd = 'ln -s %s %s' % (SIMRNA_DATA_PATH, os.getcwd())
    print('getting SimRNA data folder in cwd ...', cmd)
    logger.info(cmd)
    os.system(cmd)


def extract(template, trafl, number_of_structures=''):
    """Run SimRNA_trafl2pdb to extract all full atom structures in the trajectory."""
    cmd = 'SimRNA_trafl2pdbs %s %s AA :%s' % (template, trafl, number_of_structures)
    # cmd = 'SimRNA_trafl2pdbs %s %s AA : ' % (template, trafl)  # like extract all :
    print('rna_simrna_extract::' + cmd)
    os.system(cmd)


def cleanup(trafl):
    """Create _<trafl> with all CG structures and ss_detected and in current directory
    keep only full atom structures.

    e.g ::

            _1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low

            [mm] sim_a99 ls _1* | head
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000001.pdb
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000001.ss_detected
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000002.pdb
            1db3ee42-d2b1-4c6e-b81c-92f972e03310_ALL_100low-000002.ss_detected

    Now, in the folder there will be only AA.pdb
    find can be used to remove files, if rm can't handle, e.g. with 10k files.
    """
    traflfn = trafl.replace('.trafl', '')
    try:
        os.mkdir('_%s' % traflfn)
    except OSError:
        pass

    # move AA to a folder, .e.g _ade_pk-35b2a2c1_ALL_top500
    cmd = 'mv -v %s*AA.pdb _%s' % (traflfn, traflfn)
    print(cmd)
    os.system(cmd)

    # mv vs trash them
    if True:
        # os.system('mv %s*.pdb _%s' % (traflfn, traflfn))
        os.system('rm %s*.pdb' % traflfn)
        # os.system('find . -name "*%s*.pdb" -print0 | xargs -0 rm ' % traflfn)
        # os.system('mv %s*.ss_detected _%s' % (traflfn, traflfn))
        os.system('rm %s*.ss_detected' % traflfn)
        # os.system('find . -name "*%s*ss_detected" -print0 | xargs -0 rm ' % traflfn)
        pass

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    get_data()
    extract(args.template, args.trafl, args.number_of_structures)
    if args.cleanup:
        cleanup(os.path.basename(args.trafl))
