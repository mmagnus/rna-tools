#!/usr/bin/env python

"""
A tool to calc inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC between two structures.

ClaRNA_play required!
https://gitlab.genesilico.pl/RNA/ClaRNA_play (internal GS gitlab server)

"""
import progressbar
import argparse
import sys
import os
import subprocess
import re
import tempfile
from multiprocessing import Pool, Lock, Value, Process
number_processes = 8

from rna_pdb_tools.utils.clarna_app import clarna_app

#counter = Value(c_int)

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

    parser.add_argument('-o',"--out_fn",
                         dest="out_fn",
                         default='inf.csv',
                         help="out csv file")

    parser.add_argument('files', help="files", nargs='+')

    return parser

def do_job(i):
    i_cl_fn = clarna_app.clarna_run(i, args.force)
    output = clarna_app.clarna_compare(target_cl_fn,i_cl_fn)
    c = input_files.index(i)
    bar.update(100* (c/float(len(input_files))))
    return output

if __name__ == '__main__':
    print 'rna_calc_inf'
    print '-' * 80

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
    out_fn = args.out_fn

    # output
    print 'target, fn, inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC'
    f = open(out_fn, 'w')
    #t = 'target:' + os.path.basename(target_fn) + ' , rmsd_all\n'
    t = 'target,fn,inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC\n'
    f.write(t)

    bar = progressbar.ProgressBar()
    p = Pool(number_processes)
    # print do_job(input_files[0])
    p.map(do_job, input_files)

    #procs = []
    #for c, i in enumerate(input_files):#, range(len(input_files))):
    #    proc = Process(target=do_job, args=(i,))
    #    procs.append(proc)
    #    proc.start()

    #    print i
    #    global methods, c
    #    f.write(re.sub('\s+', ',', output) + '\n')
    #    bar.update(c/float(len(input_files)))
    #print 'csv was created! ', out_fn
