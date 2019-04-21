#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SimRNA cluster - a simple wrapper around clustering of SimRNA.

By default: cluster 1% and give 5 clusters

.. warning :: clustering has to be visiable in your shell by the script!"""

from __future__ import print_function
from rna_tools.tools.simrna_trajectory.simrna_trajectory import SimRNATrajectory
import argparse
import os


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('trafl', help="SimRNA trafl file")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    fn = args.trafl
    # 12.5
    # How to get the length?
    t = SimRNATrajectory()
    t.load_from_file(fn, only_first_frame=True)

    # print len(t) # get length of trajectory (# of frames)
    print('len of the first frame:', len(t[0]))  # get length of 1 frames (# of residues)

    # fraction of lowest energy frames to clustering: 0.01
    # rmsd thrs for clustering 0.1*seq_lenght which is: 12.5
    cutoff = 0.1 * len(t[0])
    cmd = "clustering " + fn + " 0.01 " + str(cutoff) + " | tee clustering.log"
    print(cmd)
    os.system(cmd)
