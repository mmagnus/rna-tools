#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""**run_rosetta** - wrapper to ROSETTA tools for RNA modeling

Based on C. Y. Cheng, F. C. Chou, and R. Das, Modeling complex RNA tertiary folds with Rosetta, 1st ed., vol. 553. Elsevier Inc., 2015.
http://www.sciencedirect.com/science/article/pii/S0076687914000524 

    ade$ rna_rosetta_cluster.py ade.out

The first number states how many processors to use for the run, while the second number is 1/6 the total number of previously generated FARNA models. If you are running on a supercomputer that only allows specific multiples of processors, use an appropriate number for the first input.

    rosetta_submit.py min_cmdline min_out 1 24

rosetta_submit.py min_cmdline min_out [1] [16] The first number states how many processors to use for each line in min_cmdline. Here, enter 1 for the first input so that the total number of processors used will be equal to the number of processors entered with the “-proc” flag in command line [12], above. The second number states the maximum time each job will be allowed to run (walltime). Start the run with the appropriate command listed by the out- put above (e.g., source qsubMPI for the Stampede cluster).

"""
CPUN = 200

import argparse
import argparse
import os
import glob
import subprocess
import math
import logging
import shutil

def get_no_structures(file):
    p = subprocess.Popen('cat ' + file + ' | grep SCORE | wc -l', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    stderr = p.stderr.read().strip()
    if stderr:
        print stderr
    return int(p.stdout.read().strip()) - 1

def min(f, take_n):
    """Min
    #-skip coord_constraints
    # + ' -params_file ' + f.replace('.out', '.params') +
    # ERROR:  -params_file not supported in rna_minimize anymore.
    # ERROR:: Exit from: src/apps/public/farna/rna_minimize.cc line: 156
    # -cst_fa_file
    #glycine_riboswitch_constraints glycine_riboswitch.params-ignore_zero_occupancy false -skip_ coord_constraints
    #cmd = 'extract_pdbs.default.linuxgccrelease -in::file::silent cluster.out' """

    cmd = "parallel_min_setup.py -silent " + f + " -tag " + f.replace('.out', '') + '_min  -proc ' + str(CPUN) + ' -nstruct ' + str(take_n) +  ' -out_folder min_out -out_script MINIMIZE "' + ' -ignore_zero_occupancy false "'
    print cmd
    os.system(cmd)
    os.system('chmod +x MINIMIZE && ./MINIMIZE')

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='ade.out')
    return parser

def run():
    """Pipline for modeling RNA"""
    args = get_parser().parse_args()
    ns = get_no_structures(args.file)
    print '# structures:', ns
    take_n = int(ns * 0.16) # 1/6
    min(args.file, take_n)
    
    #cluster_loop(ns)
    #cluster()
    #print get_no_structures('cluster.out')
    #if args.file:
    if 0:
        n = get_no_structures(args.file)
        nc = int(math.ceil(n * 0.005)) # nc no for clustring
        print n, nc
        get_selected(args.file, nc)
        ns = get_no_structures('selected.out')
        print ns

#main
if __name__ == '__main__':
    run()
