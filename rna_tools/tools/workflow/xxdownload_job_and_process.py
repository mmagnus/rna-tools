#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Caution: the file was changed, it does not download now from simrnaweb.
"""
from __future__ import print_function
import os

import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('--dir', action="store_false")
    parser.add_argument("rna", help="", default="evox.py -c ade")
    parser.add_argument("simrnawebid", help="", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    dir = 'simrna_' + args.simrnawebid
    os.system('mkdir  ' + dir)
    os.chdir(dir)
    #os.system('rna_simrnaweb_download_job.py --prefix tar --web-models ' + args.simrnawebid)
    #os.system('rna_simrnaweb_download_job.py  -n 200 --web-models -w ' + args.simrnawebid)
    os.system('mkdir reps')
    #os.system('mkdir reps_ns')
    os.system('cp *.pdb reps')
    #os.system('cp *.pdb reps_ns')
    os.system('evox.py -c ' + args.rna)
