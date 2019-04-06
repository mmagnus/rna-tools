#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""pefx.py - per folder eXectue cmd.

Usage::

    pefx.py 'ls | grep '.pdb\$' | wc -l'

"""
from __future__ import print_function

import argparse
import os
import glob
import sys
import re


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-d', "--dryrun",
                        action="store_true", help="dry run", default=False)

    parser.add_argument('-p', '--path', help="", default='')
    parser.add_argument('-f', '--folder-only', help="",  action="store_true")
    parser.add_argument('-c', '--case', help="only one case, for test")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('cmd')
    return parser


def sort_nicely(l):
   """ Sort the given list in the way that humans expect.
   http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
   """
   convert = lambda text: int(text) if text.isdigit() else text
   alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
   l.sort( key=alphanum_key )
   return l


def exe(cmd, dryrun, verbose):
    if verbose: print(cmd)
    if not dryrun: os.system(cmd)


def main(dryrun, path, case, cmd, folder_only):
    if path:
        os.chdir(path)
    root = os.getcwd()
    cases = sort_nicely(glob.glob('*'))

    for c in cases:
        if folder_only:
            if not os.path.isdir(c):
                continue
        #print('Case: ', c)
        #print(c, end='', flush=True)
        sys.stdout.write(c)
        sys.stdout.flush()
        # mode only for a specific case
        if case: # only if this is used
            if c != case:
                print('!!! skip ' + c + '!!!')
                continue

        # if not break ed, use this
        if os.path.isdir(c):
            # print ('Inside %s' % c)
            os.chdir(c) # go inside a folder
            # cmd
            #subcases = glob.glob('*.fa')
            #for sc in subcases:
            exe(cmd, dryrun, verbose=False)
            os.chdir(root)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    main(args.dryrun, args.path, args.case, args.cmd, args.folder_only)
