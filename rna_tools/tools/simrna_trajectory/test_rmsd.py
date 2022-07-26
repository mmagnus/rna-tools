#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse

from simrna_trajectory import SimRNATrajectory

t = SimRNATrajectory()
#t.load_from_file("Sreyoshi_seq-221ca19d_ALL_thrs1.90A_clust01.trafl")#hi_seq-221ca19d_ALL_thrs1.90A_clust01.trafl")##subset.trafl")
t.load_from_file("subset.trafl")
for y in t:
    for x in t:
        rmsd = t[0].rmsd_to(f)
        print(rmsd, )
