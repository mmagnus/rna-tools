#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
/home/magnus/rna-evo-malibu/ade # run here not up
"""
from __future__ import print_function

import argparse
import os
import glob


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-d', "--dryrun",
                        action="store_true", help="dry run", default=False)

    parser.add_argument('-p', '--path', help="", default='')
    parser.add_argument('-c', '--case', help="only one case, for test")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")

    return parser

def exe(cmd, dryrun):
    print(cmd)
    if not dryrun: os.system(cmd)

def main(dryrun, path, case):
    if path:
        os.chdir(path)
    root = os.getcwd()
    cases = glob.glob('*')
    for c in cases:
        print('Case: ', c)
        # mode only for a specific case
        if case: # only if this is used
            if c != case:
                print('!!! skip ' + c + '!!!')
                continue

        # if not break ed, use this
        if os.path.isdir(c):
            print ('Inside %s' % c)
            os.chdir(c) # go inside a folder
            # cmd
            #subcases = glob.glob('*.fa')
            #for sc in subcases:
            for i in [1000]: #
                exe('rm *min.out.*.pdb', dryrun)
                exe('mkdir %s_top%i' % (c, i), dryrun)
                exe('extract_lowscore_decoys.py *min.out %i' % (i), dryrun)
                exe('mv -v *min.out.*.pdb %s_top%i' % (c, i), dryrun)
            os.chdir(root)

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    main(args.dryrun, args.path, args.case)
