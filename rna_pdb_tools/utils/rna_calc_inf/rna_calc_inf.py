#!/usr/bin/env python

"""A tool to calc inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC between two structures.

ClaRNA_play required!
https://gitlab.genesilico.pl/RNA/ClaRNA_play (internal GS gitlab server)

import progressbar (in version 2) is required!
"""
import progressbar
import argparse
import sys
import os
import subprocess
import re
import tempfile
import csv
import shutil
from multiprocessing import Pool, Lock, Value, Process
number_processes = 7

from rna_pdb_tools.utils.clarna_app import clarna_app

def get_parser():
    parser =  argparse.ArgumentParser()#usage="%prog [<options>] <pdb files (test_data/*)>")

    parser.add_argument('-t',"--target_fn",
                           dest="target_fn",
                         default='',
                         help="pdb file")

    parser.add_argument('-s',"--ss",
                         dest="ss",
                         default='',
                         help="A:(([[))]]")

    parser.add_argument('-f',"--force",
                         dest="force",
                         action="store_true",
                         help="force to run ClaRNA")


    parser.add_argument('-v',"--verbose",
                         dest="verbose",
                         action="store_true",
                         help="be verbose")


    parser.add_argument('-o',"--out_fn",
                         dest="out_fn",
                         default='inf.csv',
                         help="out csv file")

    parser.add_argument('files', help="files", nargs='+')

    return parser

# Prepare the lock and the counter for MP
from ctypes import c_int
lock = Lock()
counter = Value(c_int)
DEBUG = False
def do_job(i):
    """Run ClaRNA & Compare, add 1 to the counter, write output 
    to csv file (keeping it locked)"""
    # run clarna & compare
    i_cl_fn = clarna_app.clarna_run(i, args.force)
    output = clarna_app.clarna_compare(target_cl_fn,i_cl_fn, DEBUG)
    if args.verbose:
        print output
    # counter and bar
    global counter
    counter.value += 1
    bar.update(counter.value)
    
    # write csv
    lock.acquire()
    csv_writer.writerow(output.split())
    csv_file.flush()
    lock.release()


#main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if len(sys.argv) == 1:
        print parser.print_help()
        sys.exit(1)

    input_files = args.files
    target_fn = args.target_fn
    ss = args.ss
    if ss:
        # generate target_fn
        target_cl_fn = clarna_app.get_ClaRNA_output_from_dot_bracket(ss, temp=False)
    else:
        target_cl_fn = clarna_app.clarna_run(target_fn, args.force)    

    # keep target save, don't overwrite it when force and
    # target is in the folder that you are running ClaRNA on
    # /tmp/tmp2nmeVB/1i6uD_M1.pdb.outCR
    d = tempfile.mkdtemp()
    tmp_target_cl_fn = d + os.sep + os.path.basename(target_cl_fn)
    shutil.copyfile(target_cl_fn, tmp_target_cl_fn)
    target_cl_fn = tmp_target_cl_fn
    ##
    
    out_fn = args.out_fn

    # Open output file
    csv_file = open(out_fn, 'w')
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerow('target,fn,inf_all,inf_stack,inf_WC,inf_nWC,SNS_WC,PPV_WC,SNS_nWC,PPV_nWC'.split(','))
    csv_file.flush()

    # Init bar and to the job
    bar = progressbar.ProgressBar(max_value=len(input_files))
    bar.update(0)
    if 1:
        p = Pool(number_processes)
        # print do_job(input_files[0])
        p.map(do_job, input_files)
    else:
        for c, i in enumerate(input_files):#, range(len(input_files))):
            do_job(i)    
    print 'csv was created! ', out_fn
