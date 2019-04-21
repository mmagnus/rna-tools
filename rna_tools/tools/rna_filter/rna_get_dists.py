#!/usr/bin/env python

"""rna_filter.py - calculate distances based on given restrants on PDB files or SimRNA trajectories.

The format of restraints::

    (d:A1-A2 < 10.0 1) = if distance between A1 and A2 lower than 10.0, score it with 1

Usage::

    $ python rna_filter.py -r test_data/restraints.txt -s test_data/CG.pdb
     d:A1-A2 10.0 measured: 6.58677550096 [x]
    test_data/CG.pdb 1.0 1 out of 1

    # $ python rna_filter.py -r test_data/restraints.txt -t test_data/CG.trafl
    (d:A1-A2 <  10.0  1)|(d:A2-A1 <= 10 1)
     restraints [('A1', 'A2', '<', '10.0', '1'), ('A2', 'A1', '<=', '10', '1')]

    Frame #1 e:1252.26
      mb for A1 [ 54.729   28.9375  41.421 ]
      mb for A2 [ 55.3425  35.3605  42.7455]
       d:A1-A2 6.58677550096
      mb for A2 [ 55.3425  35.3605  42.7455]
      mb for A1 [ 54.729   28.9375  41.421 ]
       d:A2-A1 6.58677550096
    # this ^ is off right now
"""

from __future__ import print_function
import logging
from rna_tools.rna_tools_logging import logger
from rna_tools.tools.rna_calc_rmsd.lib.rmsd.calculate_rmsd import get_coordinates
from rna_tools.tools.extra_functions.select_fragment import select_pdb_fragment_pymol_style, select_pdb_fragment
from rna_tools.tools.simrna_trajectory.simrna_trajectory import SimRNATrajectory

import argparse
import re
import numpy as np

logger.setLevel(logging.DEBUG)
logger.propagate = False


class RNAFilterErrorInRestraints(Exception):
    pass


def parse_logic(restraints_fn, verbose):
    """Parse logic of restraints.

    Args:
       restraints_nf(string): path to a file with restraints in the rigth format (see below)
       verbose (bool)       : be verbose?

    Format::

        # ignore comments
        (d:A1-A2 <  10.0  1)|(d:A2-A1 <= 10 1)

    Returns:
       list: parse restraints into a list of lists, e.g. [('A9', 'A41', '10.0', '1'), ('A10', 'A16', '10', '1')]

    """

    txt = ''
    with open(restraints_fn) as f:
        for l in f:
            if not l.startswith('#'):
                txt += l.strip()
    if verbose:
        logger.info(txt)
    restraints = re.findall(
        '\(d:(?P<start>.+?)-(?P<end>.+?)\s*(?P<operator>\>\=|\=|\<|\<\=)\s*(?P<distance>[\d\.]+)\s+(?P<weight>.+?)\)', txt)
    return restraints


def parse_logic_newlines(restraints_fn, offset=0, verbose=False):
    """Parse logic of restraints.

    Args:
       restraints_nf(string): path to a file with restraints in the rigth format (see below)
       verbose (bool)       : be verbose?

    Format::

        # ignore comments
        d:Y23-Y69 < 25.0
        d:Y22-Y69 < 25.0
        # d:<chain><resi_A>-<resi_B> <operator> <distance> <weight>; each restraints in a new line

    Raises:
       __main__.RNAFilterErrorInRestraints: Please check the format of your restraints!

    Returns:
       list: parse restraints into a list of lists, e.g. [('A9', 'A41', '10.0', '1'), ('A10', 'A16', '10', '1')]

    """
    restraints = []
    with open(restraints_fn) as f:
        for l in f:
            if l.strip():
                if not l.startswith('#'):
                    if verbose:
                        logger.info(l)
                    restraint = re.findall(
                        'd:(?P<start>.+?)-(?P<end>.+?)\s*(?P<operator>\>\=|\=|\<|\<\=)\s*(?P<distance>[\d\.]+)\s+(?P<weight>.+?)', l)
                    if restraint:
                        # without [0] it is restraints [[('Y23', 'Y69', '<', '25.0', '1')], [('Y22', 'Y69', '<', '25.0', '1')]]
                        # why? to convert 'Y23', 'Y69', '<', '25.0', '1' -> 'Y23', 'Y69', '<', 25.0, 1
                        start = restraint[0][0][0] + str(int(restraint[0][0][1:]) + offset)
                        end = restraint[0][1][0] + str(int(restraint[0][1][1:]) + offset)
                        restraints.append([start, end, restraint[0][1], restraint[0][2],
                                           float(restraint[0][3]), float(restraint[0][4])])

    if len(restraints) == 0:
        raise RNAFilterErrorInRestraints('Please check the format of your restraints!')
    return restraints  # [('A9', 'A41', '10.0', '1'), ('A10', 'A16', '10', '1')]


def get_distance(a, b):
    diff = a - b
    return np.sqrt(np.dot(diff, diff))


def parse_pdb(pdb_fn, selection):
    """
{'A9': {'OP1': array([ 53.031,  21.908,  40.226]), 'C6': array([ 54.594,  27.595,  41.069]), 'OP2': array([ 52.811,  24.217,  39.125]), 'N4': array([ 53.925,  30.861,  39.743]), "C1'": array([ 55.611,  26.965,  43.258]), "C3'": array([ 53.904,  25.437,  43.809]), "O5'": array([ 53.796,  24.036,  41.353]), 'C5': array([ 54.171,  28.532,  40.195]), "O4'": array([ 55.841,  25.746,  42.605]), "C5'": array([ 54.814,  23.605,  42.274]), 'P': array(
    [ 53.57 ,  23.268,  39.971]), "C4'": array([ 55.119,  24.697,  43.283]), "C2'": array([ 54.563,  26.706,  44.341]), 'N1': array([ 55.145,  27.966,  42.27 ]), "O2'": array([ 55.208,  26.577,  45.588]), 'N3': array([ 54.831,  30.285,  41.747]), 'O2': array([ 55.76 ,  29.587,  43.719]), 'C2': array([ 55.258,  29.321,  42.618]), "O3'": array([ 53.272,  24.698,  44.789]), 'C4': array([ 54.313,  29.909,  40.572])}}
    """
    V = {}
    with open(pdb_fn) as f:
        for line in f:
            if line.startswith("ATOM"):
                curr_chain_id = line[21]
                curr_resi = int(line[22: 26])
                curr_atom_name = line[12: 16].strip()
                if selection:
                    if curr_chain_id in selection:
                        if curr_resi in selection[curr_chain_id]:
                            x = line[30: 38]
                            y = line[38: 46]
                            z = line[46: 54]
                            # V.append(np.asarray([x,y,z],dtype=float))
                            if curr_chain_id + str(curr_resi) in V:
                                V[curr_chain_id +
                                    str(curr_resi)][curr_atom_name] = np.asarray([x, y, z], dtype=float)
                            else:
                                V[curr_chain_id + str(curr_resi)] = {}
                                V[curr_chain_id +
                                    str(curr_resi)][curr_atom_name] = np.asarray([x, y, z], dtype=float)
    return V


def check_condition(condition, wight):
    """return True/False, score"""
    pass


def get_residues(pdb_fn, restraints, verbose):
    residues = set()
    for h in restraints:
        a = h[0]
        b = h[1]
        a = a[0] + ':' + a[1:]
        residues.add(a)  # A19
        b = b[0] + ':' + b[1:]
        residues.add(b)
    # set(['A:41', 'A:9', 'A:10', 'A:16'])

    selection = ','.join(residues)
    selection_parsed = select_pdb_fragment(selection, separator=",", splitting="[,:;]")

    residues = parse_pdb(pdb_fn, selection_parsed)

    # get mb
    for r in residues:
        if 'N9' in residues[r]:  # A,G
            residues[r]['mb'] = residues[r]['N9'] - ((residues[r]['N9'] - residues[r]['C6']) / 2)
        else:  # A,G
            residues[r]['mb'] = residues[r]['N1'] - ((residues[r]['N1'] - residues[r]['C4']) / 2)
    for r in residues:
        if verbose:
            logger.info(' '.join(['mb for ', str(r), str(residues[r]['mb'])]))
    return residues


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-r', "--restraints_fn",
                        dest="restraints_fn",
                        required=True,
                        help="""restraints_fn:
Format:
(d:A9-A41 <  10.0  1)|(d:A41-A9 <= 10 1)
""")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")

    parser.add_argument('-s', dest="structures", help='structures',
                        nargs='+')  # , type=string)
    parser.add_argument(
        '--offset', help='use offset to adjust your restraints to numbering in PDB files, ade (1y26)'
        'pdb starts with 13, so offset is -12)', default=0, type=int)
    parser.add_argument('-t', dest="trajectory", help="SimRNA trajectory")
    return parser


def calc_dists_for_pdbs(pdb_files, pairs, verbose):
    """
    """
    # h = ('A1', 'A2', '<', '10.0', '1')
    for pdb_fn in pdb_files:
        # logger.info(pdb_fn)
        score = 0
        residues = get_residues(pdb_fn, restraints, verbose)
        good_dists = 0
        for h in pairs:
            dist = get_distance(residues[h[0]]['mb'], residues[h[1]]['mb'])
            # change distance
            ok = '[ ]'
            if dist < h[4]:
                score += h[5]
                ok = '[x]'
                good_dists += 1
            print(' '.join([' d:' + h[0] + '-' + h[1] + ' ' + str(h[4]), 'measured:', str(dist), ok]))
        print(pdb_fn, score / float(len(restraints)), good_dists, 'out of', len(restraints))


# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    pairs = eval(args.pairs)
    calc_dists_for_pdbs(args.structures, pairs, args.verbose)
