#!/usr/bin/env python

"""rna_align_find_core.py
"""

import rna_tools.tools.rna_alignment.rna_alignment as ra
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="fasta seq")
    return parser

if __name__ == '__main__':
    # for debugging #
    #fn = 'test_data/RF00167.stockholm.sto'
    #ids = ['AL591975.1/251136-251218','CP000721.1/2204691-2204778']#ACCL02000010.1/116901-116991']
    parser = get_parser()
    args = parser.parse_args()
    fn = args.file

    a = ra.RNAalignment(fn)
    print(a.find_core())
