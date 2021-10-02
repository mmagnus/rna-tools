#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from rna_tools.rna_tools_lib import RNAStructure
from rna_tools.tools.md.rnakb_utils import make_rna_gromacs_ready
import argparse

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

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        r = RNAStructure(f)
        for i in r.get_all_chain_ids():
            chain = r.get_chain(i)
            nchain = make_rna_gromacs_ready(chain)
            print(nchain)
