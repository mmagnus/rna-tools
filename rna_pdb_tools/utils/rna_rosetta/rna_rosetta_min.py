#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_rosetta_min.py - a script to do minimization

The script takes the number of structures and the analyzed silence file and does the maths.

Job names will be as your silent file preceding with ~, .e.g ``~tha``.

http://www.sciencedirect.com/science/article/pii/S0076687914000524 ::

    ade$ rna_rosetta_cluster.py ade.out

The first number states how many processors to use for the run, while the second number is 1/6 the total number of previously generated FARNA models. If you are running on a supercomputer that only allows specific multiples of processors, use an appropriate number for the first input.::

    rosetta_submit.py min_cmdline min_out 1 24

rosetta_submit.py min_cmdline min_out [1] [16] The first number states how many processors to use for each line in min_cmdline. Here, enter 1 for the first input so that the total number of processors used will be equal to the number of processors entered with the "-proc" flag in command line [12], above. The second number states the maximum time each job will be allowed to run (walltime). Start the run with the appropriate command listed by the out- put above (e.g., source qsubMPI for the Stampede cluster).

E.g. for 20k silet file, 1/6 will be minimized = 3.3k::

    parallel_min_setup.py -silent rp21cr62.out -tag rp21cr62_min  -proc 200 -nstruct 3200 -out_folder mo -out_script MINIMIZE " -ignore_zero_occupancy false "
    rosetta_submit.py MINIMIZE mo 1 100 m

    [peyote2] rp21 easy_cat.py mo
    Catting into:  rp21_min.out ... from 200 primary files. Found 3200  decoys.

    # on 200 cpus it took around ~30min
"""
from __future__ import print_function
import logging

logging.basicConfig(filename='rna_rosetta_min.log', level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s")

import argparse
import os
import glob
import subprocess
import math
import logging
import shutil
import re


def get_no_structures(file):
    """Get a number of structures in a silent file"""
    p = subprocess.Popen('cat ' + file + ' | grep SCORE | wc -l', shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    stderr = p.stderr.read().strip()
    if stderr:
        print(stderr)
    return int(p.stdout.read().strip()) - 1


def min(silent_file, take_n, cpus, go):
    """Run parallel_min_setup (to MINIMIZE file), rosetta_submit.py, and qsubMINI.

    Fix on the way, qsub files::

        -out:file:silent mo/0/mo/123/tha_min.out -> -out:file:silent mo/123/tha_min.out

    I don't know why mo/0/ is there. I might be because of my changes in rosetta_submit.py (?). """

    #-skip coord_constraints
    # + ' -params_file ' + f.replace('.out', '.params') +
    # ERROR:  -params_file not supported in rna_minimize anymore.
    # ERROR:: Exit from: src/apps/public/farna/rna_minimize.cc line: 156
    # -cst_fa_file
    # glycine_riboswitch_constraints glycine_riboswitch.params-ignore_zero_occupancy false -skip_ coord_constraints
    #cmd = 'extract_pdbs.default.linuxgccrelease -in::file::silent cluster.out'

    logging.info('cpun %i' % cpus)

    # parallel_min_setup
    cmd = "parallel_min_setup.py -silent " + silent_file + " -tag " + silent_file.replace('.out', '') + '_min  -proc ' + str(
        cpus) + ' -nstruct ' + str(take_n) + ' -out_folder mo -out_script MINIMIZE "' + ' -ignore_zero_occupancy false "'
    print(cmd)
    logging.info(cmd)
    os.system(cmd)

    # rosetta_submit
    cmd = "rosetta_submit.py MINIMIZE mo 1 100 m"
    print(cmd)
    logging.info(cmd)
    os.system(cmd)

    # fixing qsub file
    d = "qsub_files"

    for c, qf in enumerate(os.listdir(d), 0):
        qf = d + os.sep + qf
        f = open(qf)
        txt = f.read()
        #x = re.search('-out:file:silent (?P<extra>mo\/\d+)\/min', txt )
        txt = txt.replace('-out:file:silent mo/0/', '-out:file:silent ')

        # fix naming of my jobs
        txt = txt.replace('#$ -N mm0', '#$ -N ~' +
                          silent_file.replace('.out', '').strip() + '_' + str(c))  # to get ~rp06bx571_186

        f.close()

        # write to file
        open(qf, 'w').write(txt)

    if go:
        # run qsubMINI
        logging.info('qsubMINI')
        os.system('chmod +x ./qsubMINI')
        os.system('./qsubMINI')


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help='ade.out')
    parser.add_argument('-g', '--go', action='store_true')
    parser.add_argument('-c', '--cpus', help='default: 200', default=200)
    return parser


def run():
    args = get_parser().parse_args()
    ns = get_no_structures(args.file)
    print('# structures:', ns)
    take_n = int(ns * 0.16)  # 1/6
    min(args.file, take_n, int(args.cpus), args.go)


# main
if __name__ == '__main__':
    run()
