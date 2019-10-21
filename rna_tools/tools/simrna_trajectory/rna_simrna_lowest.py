#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_simrna_lowest.py - get lowest energy frames out of a SimRNA trajectory file

This code uses heavily the SimRNATrajectory class. Be default 100 lowest energy frames is exported.
"""

from __future__ import print_function
from rna_tools.tools.simrna_trajectory.simrna_trajectory import SimRNATrajectory
import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--nstruc', help="SimRNA trafl file", type=int, default=100)
    parser.add_argument('trafl', help="SimRNA trafl file")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    s = SimRNATrajectory()
    fn = args.trafl

    s.load_from_file(fn, top_level=True)

    # for f in s.frames:
    #    print f

    sorted_frames = s.sort(inplace=False)

    for c, f in enumerate(sorted_frames[:args.nstruc]):
        print(c + 1, f)
        # print f.header
        # print f.coords

    s2 = SimRNATrajectory()
    s2.load_from_list(sorted_frames[:args.nstruc])
    s2.plot_energy(fn.replace('.trafl', '.png'))
    s2.save(fn.replace('.trafl', '') + '_top' + str(args.nstruc) + '.trafl')
