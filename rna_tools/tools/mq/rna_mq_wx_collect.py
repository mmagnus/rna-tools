#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
/home/magnus/rna-evo-malibu/ade # run here not up

at VM
[7]  + 15507 suspended  rna_pefx.py 'rna_mq_collect.py -t FARFAR2_hires -g FARFAR2*csv struc/*'


"""
from __future__ import print_function

import argparse
import os
import glob

import re
def sort(l):
    """ Sort the given list in the way that humans expect.
    http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
    """
    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', key)]
    l.sort(key=alphanum_key)
    return l


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-d', "--dryrun",
                        action="store_true", help="dry run", default=False)
    parser.add_argument('-p', '--path', help="", default='')
    parser.add_argument('--dataset', help="only one case, for test", default="*")
    parser.add_argument('-e', '--exe', help="only one case, for test",
                        default="rna_mq_collect.py -t FARFAR2_hires struc/*.pdb -m 4 -v -f -o FARFAR2_hires.csv")# -o FARFAR2")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    return parser

def exe(cmd, dryrun=False):
    print(cmd)
    if not dryrun: os.system(cmd)

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if args.path:
        os.chdir(args.path)
    root = os.getcwd()
    dataset = sort(glob.glob(args.dataset))# '*')
    
    for c in dataset:
        print('dataset ', c)
        if c.startswith('_'):
            continue
        # mode only for a specific case
        #if args.case: # only if this is used
        #     if c != args.case:
        #        print('!!! skip ' + c + '!!!')
        #        continue
        # if not break ed, use this
        if os.path.isdir(c):
            os.chdir(root + '/' + c)
            subcases = sort(glob.glob('*'))
            for sc in subcases:
                os.chdir(sc) # go inside a folder
                print ('- case %s' % sc)
                if sc.startswith('_'):
                    continue
                if not args.dryrun:
                    exe(args.exe)
                os.chdir('..')
            ## for i in [1000]: #
            ##     exe('rm *min.out.*.pdb', dryrun)
            ##     exe('mkdir %s_top%i' % (c, i), dryrun)
            ##     exe('extract_lowscore_decoys.py *min.out %i' % (i), dryrun)
            ##     exe('mv -v *min.out.*.pdb %s_top%i' % (c, i), dryrun)
        os.chdir(root)

