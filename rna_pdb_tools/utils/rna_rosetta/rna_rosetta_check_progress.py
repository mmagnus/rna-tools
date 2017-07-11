#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_rosetta_check_progress.py

Usage::

    [peyote2] rosetta_jobs rna_rosetta_check_progress.py .
             jobs  #curr  #todo  #decoys done
    0  ./rp17s223    200      0      407  [ ]
    1   ./rp17hcf      0      0        0  [ ]
    #curr  232 #todo  0

TODO: Check progress and kill!
"""

import commands
import os
import pandas as pd
import glob
import sys
import argparse

try:
    from rna_pdb_tools.rpt_config import RNA_ROSETTA_RUN_ROOT_DIR_MODELING, EASY_CAT_PATH
except:
    print ('Set up rna_rosetta_run_root_dir_for_modeling in rpt_config_local.py')

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', help='', default=RNA_ROSETTA_RUN_ROOT_DIR_MODELING)
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args()
    v = args.verbose
    
    #jobs = [x.replace('.fa','') for x in glob.glob( '*.fa')] # sys.argv[1] +
     # sys.argv[1] +
    #print jobs

    jobs = [x for x in glob.glob(sys.argv[1] + '/*')]
    if v: print(jobs)
    
    curr_path = os.getcwd()
    limit = 2000 # how many decouys you want

    d = {}
    d['#decoys'] = []
    d['jobs'] = []
    d['#curr'] = []
    d['#todo'] = []
    d['done'] = []


    for j in jobs:
        if v: print(j)
        d['jobs'].append(j)
        try:
            os.chdir(j)
        except:
            pass # OSError: [Errno 20] Not a directory:

        dirs = []
        # Uff.. this is crazy.
        for c in ['_', '__', '___', 'x', 'y', 'z', 'X', 'Y', 'Z', 'o', 'out'] + ['r' + str(i) for i in range(0, 1000)]:
            if os.path.exists(c):
                dirs.append(c)
        # print ' ', j, '', dirs
        cmd = EASY_CAT_PATH + ' ' + ' '.join(dirs)
        if v: print(cmd)
        out = commands.getoutput(cmd)
        if out.strip():
            if len(j) <= 5:
                if v: print out + '\t\t',
            else:
                if v: print out + '\t',
            no_decoys = int(out.split()[-2])
        else:
            no_decoys = 0       

        if v: print '#', no_decoys

        # check current running, only 6 char are taken from dir name
        # /home/magnus/rosetta_jobs/rp17s223 -> rp17s2199, rp17s2185
        cmd = 'qstat | grep ' + os.path.basename(j)[:6] + ' | grep "  r  " | wc -l '
        if v: print(cmd)
        out = commands.getoutput(cmd) # aacy97r <- r
        if v:print '@cluster #curr ', out,
        curr = int(out)
        d['#curr'].append(curr)

        # check todo
        out = commands.getoutput('qstat | grep ' + os.path.basename(j)[:6] + ' | grep "  qw  " | wc -l ')
        if v:print '#todo ', out,
        todo = int(out)
        d['#todo'].append(todo)

        d['#decoys'].append(no_decoys)

        if no_decoys > limit:
            if v:print '# @cluster:', out, ' < ............... OK'
            d['done'].append('[✔]')
            ## killing
            #cmd = 'kill '
            #print cmd
            #os.system('
        else:
            d['done'].append('[ ]')
            if v:print '# @cluster:', out
            ##cmd = 'kill '
            #print cmd

        os.chdir(curr_path)

    df = pd.DataFrame(d, columns=['jobs', '#curr', '#todo', '#decoys', 'done'])
    print(df)

    out = commands.getoutput('qstat | grep magnus  | grep "  r  " | wc -l ')
    print '#curr ', out,
    out = commands.getoutput('qstat | grep magnus | grep "  qw  " | wc -l ')
    print '#todo ', out,

