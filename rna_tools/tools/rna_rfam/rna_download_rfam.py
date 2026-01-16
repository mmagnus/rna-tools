#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_download_rfam.py - download all cm and seeds

For the list grep seed file::

    (base) rnahub@rnahub:~/db/rfam$ grep '#=GF AC' Rfam.seed  | wc -l
    #=GF AC   RF04297
    #=GF AC   RF04298
    #=GF AC   RF04299
    #=GF AC   RF04300

or cm::

    (base) rnahub@rnahub:~/db/rfam$ grep ACC Rfam.cm
    ACC      RF04235
    ACC   RF04235
    ACC      RF04236
    ACC   RF04236


"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import os

def process_rfam(family_id, job_path='rfam'):
        cmd = f"wget https://rfam.org/family/{family_id}/alignment/stockholm -O {job_path}/{family_id}.seed.sto"
        print(cmd)
        os.system(cmd)
        cmd = f"wget https://rfam.org/family/{family_id}/alignment/fastau -O {job_path}/{family_id}.seed.fa"
        print(cmd)
        os.system(cmd)
        cmd = f"wget https://rfam.org/family/{family_id}/cm -O {job_path}/{family_id}.cm"
        print(cmd)
        os.system(cmd)

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    for l in open(args.file):
        rfam_id = l.split()[2]
        print(rfam_id)
        process_rfam(rfam_id)
