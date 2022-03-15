#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_rosetta_extract_lowscore_decoys.py - a simple wrapper to extract_lowscore_decoys.py

To be used in Jupter notebooks and other scripts.
"""
from __future__ import print_function

import argparse
import subprocess
import os
import logging
import time

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("nstruc", help="# of low score structures to obtained", default=100)
    parser.add_argument("file", help="silent file")
    return parser


def rosetta_extract_lowscore_decoys(file, nstruc, log=True):
    cmd = 'extract_lowscore_decoys.py ' + file + ' ' + nstruc
    os.system(cmd)

    if log:
        logging.basicConfig(filename='rna_rosetta_extract_lowscore_decoys.log',level=logging.INFO)
        #logging.basicConfig(level=logging.INFO)
        logging.info(args)
        logging.info(time.strftime("%Y-%m-%d %H:%M"))
        logging.info(args)

    ## p = subprocess.Popen(cmd,
    ##                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    ## stderr = p.stderr.read().strip()
    ## if stderr:
    ##     print(stderr)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    rosetta_extract_lowscore_decoys(args.file, args.nstruc)
