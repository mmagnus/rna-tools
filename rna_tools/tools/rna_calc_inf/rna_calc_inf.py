#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A tool to calc inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC between two structures.

Mind, that ClaRNA is pretty slow, it takes even a few seconds to analyze a structure,
so for, say, 1000 models you need a few hours.

How to make it faster? First, you can use ``--number_of_threads`` to specify the number of cores used for multiprocessing.

Second, the procedure implemented in here is composed of two steps, first for each structure ClaRNA is used to generate an output with contacts, then these files are used for comparisons. So, if you want to re-run your analysis, you don't have to run re-run ClaRNA itself. Thus, be default ClaRNA is not executed if <model>.outCR is found next to the analyzed files.  To change this behavior force (``--force``) rna_cal_inf.py to re-run ClaRNA.

"""
from __future__ import print_function

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
from rna_tools.tools.clarna_app import clarna_app
#from rna_tools.opt.BasicAssessMetrics.BasicAssessMetrics import InteractionNetworkFidelity

import pandas as pd
pd.set_option('display.max_rows', 1000)


def get_parser():
    parser =  argparse.ArgumentParser()#usage="%prog [<options>] <pdb files (test_data/*)>")
    parser.add_argument('-t',"--target_fn",
                           dest="target_fn",
                         default='',
                         help="pdb file")

    parser.add_argument('-m',"--number_of_threads",
                         dest="nt",
                         default=3,
                         help="number of threads used for multiprocessing, if 1 then mp is not used \
                         (useful for debugging)!")

    parser.add_argument('--ignore-files', help='files to be ingored, .e.g, \'solution\'', default='')

    parser.add_argument('-s',"--ss",
                         dest="ss",
                         default='',
                         help="A:(([[))]], works only for single chain (the chain is A by default)")

    parser.add_argument('--no-stacking',
                         action="store_true",
                         help="default: use stacking, if this option on, don't take into account stacking, \n\nWARNING/BUG: inf_all will be incorrectly calculated if stacking is off")

    parser.add_argument('--debug',
                         action="store_true")

    parser.add_argument('-pr', '--print-results',
                         action="store_true")

    parser.add_argument('-sr', '--sort-results',
                         action="store_true")

    parser.add_argument('--method', default="clarna", help="you can use mcannotate* or clarna (right now only clarna is tested)")

    parser.add_argument("--target-selection",
                         default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")

    parser.add_argument("--model-selection",
                         default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")

    parser.add_argument('--renumber-residues', help='renumber residues from 1 to X for comparison with selection',
                        default='',
                        action="store_true")

    parser.add_argument('--dont-remove-sel-files',
                         action="store_true",
                         help="don't remove temp files created based on target|model-selectionforce")
     
    parser.add_argument('-f',"--force",
                         dest="force",
                         action="store_true",
                         help="force to run ClaRNA even if <pdb>.outCR file is there, for will be auto True when selection defined")

    parser.add_argument('-v',"--verbose",
                         dest="verbose",
                         action="store_true",
                         help="be verbose, tell me more what're doing")

    parser.add_argument('-o',"--out_fn",
                         dest="out_fn",
                         default='inf.csv',
                         help="out csv file, be default `inf.csv`")

    parser.add_argument('files', help="files, .e.g folder_with_pdbs/*pdbs", nargs='+')
    return parser

# Prepare the lock and the counter for MP
from ctypes import c_int
lock = Lock()
counter = Value(c_int)

def do_job(i, method='clarna'):
    """Run ClaRNA & Compare, add 1 to the counter, write output
    to csv file (keeping it locked)"""
    #if method == 'clarna':
        # run clarna & compare
    i_cl_fn = clarna_app.clarna_run(i, args.force,not args.no_stacking)
    output = clarna_app.clarna_compare(target_cl_fn,i_cl_fn, verbose=DEBUG)
    if args.verbose:
        print(output)
    ## else:
    ##     rmsd, DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK = InteractionNetworkFidelity(os.path.abspath(target_fn),
    ##                                                                                   '/tmp/empty-index',
    ##                                                                                   os.path.abspath(i),
    ##                                                                                   '/tmp/empty-index')
    ##     if args.debug:
    ##         print(rmsd)

    # counter and bar
    global counter
    counter.value += 1
    bar.update(counter.value)

    # write csv
    lock.acquire()
    # take only filename of target
    cells = output.split()
    cells[0] = os.path.basename(cells[0])

    csv_writer.writerow(cells)
    csv_file.flush()
    lock.release()

#main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    DEBUG = False
    if args.debug:
        DEBUG = True
        args.verbose = True
        print(args)
        
    if len(sys.argv) == 1:
        print((parser.print_help()))
        sys.exit(1)

    input_files = args.files
    ####################################
    # implement ignore files
    tmp = []
    if args.ignore_files:
        for f in input_files:
            if args.ignore_files in f:
                continue
            tmp.append(f)
        input_files = tmp

    if args.model_selection or args.target_selection:
        args.force = True
        
    if args.model_selection:
        tmp = []
        for f in input_files:
            new_f = f.replace('.pdb', '_sel.pdb')
            cmd =  "rna_pdb_toolsx.py --extract '" + args.model_selection + "' " + f + ' > ' + new_f
            os.system(cmd)
            if args.renumber_residues:
                cmd = "rna_pdb_toolsx.py --no-hr --rpr --inplace --renumber-residues " + new_f
                os.system(cmd)
            tmp.append(new_f)
        input_files = tmp

    target_fn = args.target_fn
    if args.target_selection:
        new_target_fn = target_fn.replace('.pdb', '_sel.pdb')
        cmd =  "rna_pdb_toolsx.py --extract '" + args.target_selection + "' " + target_fn + ' > ' + new_target_fn
        os.system(cmd)
        if args.renumber_residues:
            cmd = "rna_pdb_toolsx.py --no-hr --rpr --inplace --renumber-residues " + new_target_fn
            os.system(cmd)
        target_fn = new_target_fn

    ss = args.ss
    if ss:
        # generate target_fn
        ss_txt = open(ss).read().split('\n')[2]
        target_cl_fn = clarna_app.get_ClaRNA_output_from_dot_bracket(ss_txt, temp=False)
    else:
        target_cl_fn = clarna_app.clarna_run(target_fn, args.force)
        
        
        

    # keep target save, don't overwrite it when force and
    # target is in the folder that you are running ClaRNA on
    # /tmp/tmp2nmeVB/1i6uD_M1.pdb.outCR
    d = tempfile.mkdtemp()
    tmp_target_cl_fn = d + os.sep + os.path.basename(target_cl_fn)
    shutil.copyfile(target_cl_fn, tmp_target_cl_fn)
    target_cl_fn = tmp_target_cl_fn

    out_fn = args.out_fn

    if args.no_stacking:
        print('WARNING/BUG: inf_all will be incorrectly calculated if stacking is off')

    # Open output file
    csv_file = open(out_fn, 'w')
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerow('target,fn,inf_all,inf_stack,inf_WC,inf_nWC,sns_WC,ppv_WC,sns_nWC,ppv_nWC'.split(','))
    csv_file.flush()

    # Init bar and to the job
    try:
        bar = progressbar.ProgressBar(max_value=len(input_files))
        bar.update(0)
    except TypeError:
        print('Please install progressbar2 (not progressbar), e.g. pip install progressbar2')
        sys.exit(1)

    # Main meat
    number_processes = int(args.nt)

    open('/tmp/empty-index', 'a').close()  ## ugly hack

    if number_processes > 1: # multi
        p = Pool(number_processes)
        p.map(do_job, input_files) # , args.method)
    else: # single process
        for c, i in enumerate(input_files):#, range(len(input_files))):
            do_job(i, args.method)
    print('\ncsv was created! ', out_fn)
    
    # hack with pandas
    csv_file.close()

    df = pd.read_csv(out_fn)
    df = df.round(2)
    df['target'] = df['target'].str.replace('.outCR', '')
    df['fn'] = df['fn'].str.replace('.outCR', '')
    if args.sort_results:
        df = df.sort_values('inf_all', ascending=False)
    if args.print_results:
        print(df)
    df.to_csv(out_fn, sep=',', index=False)

    # remove temp files
    if not args.dont_remove_sel_files:
        if args.target_selection:
            os.remove(os.path.abspath(target_fn))
        if args.model_selection:
            for f in input_files:
                if f == target_fn:
                    continue
                os.remove(os.path.abspath(f))
