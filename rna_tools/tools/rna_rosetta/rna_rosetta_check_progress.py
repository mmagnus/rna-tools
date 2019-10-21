#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_rosetta_check_progress.py - check progress for many simulations of Rosetta

Example::

    [peyote2] rosetta_jobs rna_rosetta_check_progress.py .
             jobs  #curr  #todo  #decoys done
    0  ./rp17s223    200      0      407  [ ]
    1   ./rp17hcf      0      0        0  [ ]
    # curr  232 #todo  0

"""
from __future__ import print_function
import subprocess

import os
import pandas as pd
import glob
import sys
import argparse

try:
    from rna_tools.rna_tools_config import RNA_ROSETTA_RUN_ROOT_DIR_MODELING, RNA_ROSETTA_NSTRUC
except:
    print('Set up RNA_ROSETTA_RUN_ROOT_DIR_MODELING in rna_tools.rna_tools_config_local.py')

from rna_tools.rna_tools_config import EASY_CAT_PATH
# print('Set up EASY_CAT_PATH in rpt_config_local.py')
print(EASY_CAT_PATH)


def run_cmd(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dir', default=RNA_ROSETTA_RUN_ROOT_DIR_MODELING,
                        help="""directory with rosetta runs, define by RNA_ROSETTA_RUN_ROOT_DIR_MODELING
right now: \n""" + RNA_ROSETTA_RUN_ROOT_DIR_MODELING)
    parser.add_argument('-v', '--verbose',
                        action='store_true', help="be verbose")
    parser.add_argument('-m', '--min-only', action='store_true',
                        help="check only for mo folder")
    parser.add_argument('-s', '--select',
                        help="select for analysis only jobs with this phrase, .e.g., evoseq_", default='')

    parser.add_argument('-k', '--kill', action='store_true', help="""kill (qdel) jobs if your reach
limit (nstruc) of structure that you want, right now is %i structures""" % RNA_ROSETTA_NSTRUC)
    return parser


if __name__ == '__main__':
    args = get_parser().parse_args()
    v = args.verbose

    # jobs = [x.replace('.fa','') for x in glob.glob( '*.fa')] # sys.argv[1] +
    # sys.argv[1] +
    # print jobs

    jobs = [x for x in glob.glob(args.dir + '/*')]
    if v:
        print(jobs)

    curr_path = os.getcwd()

    d = {}
    d['#decoys'] = []
    d['jobs'] = []
    d['#curr'] = []
    d['#todo'] = []
    d['done'] = []
    d['progress'] = []

    for j in jobs:
        if args.select:
            if not args.select in j:
                continue
        if v:
            print(j)
        d['jobs'].append(j)
        try:
            os.chdir(j)
        except:
            pass  # OSError: [Errno 20] Not a directory:

        dirs = []

        if not args.min_only:
            # Uff.. this is crazy.
            for c in ['_', '__', '___', 'x', 'y', 'z', 'X', 'Y', 'Z', 'o', 'out'] + ['r' + str(i) for i in range(0, 1000)]:
                if os.path.exists(c):
                    dirs.append(c)
            # print ' ', j, '', dirs
            cmd = EASY_CAT_PATH + ' ' + ' '.join(dirs)
            if v:
                print('cmd: %s' % cmd)

            out, err = run_cmd(cmd)

            if out.strip():
                if len(j) <= 5:
                    if v:
                        print(out + '\t\t', end='')
                else:
                    if v:
                        print(out + '\t', end='')
                no_decoys = int(out.split()[-2])
            else:
                no_decoys = 0

            if v:
                print('#', no_decoys)
        else:
            # # of minimized
            RNA_ROSETTA_NSTRUC = 1600  # change expected !!!! pretty ugly
            for c in ['mo']:
                if os.path.exists(c):
                    dirs.append(c)
            # print ' ', j, '', dirs
            cmd = EASY_CAT_PATH + ' ' + ' '.join(dirs)
            if v:
                print(cmd)

            out, err = run_cmd(cmd)

            if out.strip():
                if len(j) <= 5:
                    if v:
                        print(out + '\t\t', end='')
                else:
                    if v:
                        print(out + '\t', end='')
                no_decoys = int(out.split()[-2])
            else:
                no_decoys = 0

        if v:
            print('#', no_decoys)

        # check current running, only 6 char are taken from dir name
        # /home/magnus/rosetta_jobs/rp17s223 -> rp17s2199, rp17s2185
        cmd = 'qstat | grep ' + \
            os.path.basename(j)[:6] + ' | grep "  r  " | wc -l '
        if v:
            print(cmd)

        out, err = run_cmd(cmd)

        if v:
            print('@cluster #curr ', out, end="")
        curr = int(out)
        d['#curr'].append(curr)

        # check todo
        out, err = run_cmd('qstat | grep ' + os.path.basename(j)
                           [:6] + ' | grep "  qw  " | wc -l ')
        if v:
            print('#todo ', out, end="")

        todo = int(out)
        d['#todo'].append(todo)

        d['#decoys'].append(no_decoys)

        if no_decoys >= RNA_ROSETTA_NSTRUC:
            if v:
                print('# @cluster:', out, ' < ............... OK')
            d['done'].append('[x]')
            if args.kill:
                cmd = 'qstat | grep ' + \
                    os.path.basename(j)[:6] + \
                    " | awk '{print $1}' | xargs qdel "
                print(cmd)
                os.system(cmd)
        else:
            d['done'].append('[ ]')
            if v:
                print('# @cluster:', out)
            # cmd = 'kill '
            # print cmd

        d['progress'].append(round((float(no_decoys) / 10000) * 100, 2))

        os.chdir(curr_path)

    df = pd.DataFrame(d, columns=['jobs', '#curr', '#todo', '#decoys', 'done', 'progress'])
    print(df)

    out, err = run_cmd('qstat | grep magnus  | grep "  r  " | wc -l ')
    print('#curr ', out, end=" ")
    out, err = run_cmd('qstat | grep magnus | grep "  qw  " | wc -l ')
    print('#todo ', out)
