#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_simrna_energy_plot.py - get an energy plot of SimRNA trajectory
"""

from __future__ import print_function
from collections import deque
from rna_tools.tools.simrna_trajectory.simrna_trajectory import SimRNATrajectory
import numpy as np
import gc
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    s = SimRNATrajectory()
    s.load_from_file(args.file, top_level=True)
    #if args.verbose:
    s.plot_energy('plot.png')
