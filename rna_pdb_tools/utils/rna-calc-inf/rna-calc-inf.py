#!/usr/bin/python

"""
ClaRNA_play required!
https://gitlab.genesilico.pl/RNA/ClaRNA_play (internal GS gitlab server)
"""

import optparse
import sys
import os
import subprocess
import re

def clarna_run(fn, force):
    fn_out = fn + '.outCR'
    if os.path.isfile(fn_out) and not force:
        pass
    else:
        cmd = 'clarna_run.py -ipdb ' + fn + ' > ' + fn_out
        print cmd
        os.system(cmd)
    return fn_out

def clarna_compare(target_cl_fn,i_cl_fn):
    cmd = 'clarna_compare.py -iref ' + target_cl_fn + ' -ichk ' + i_cl_fn
    o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std = o.stdout.read().strip()
    return std
    
if __name__ == '__main__':
    print 'rna-calc-inf'
    print '-' * 80
    
    optparser=optparse.OptionParser(usage="%prog [<options>] <pdb files (test_data/*)>")

    optparser.add_option('-t',"--target_fn", type="string",
                         dest="target_fn",
                         default='',
                         help="pdb file")

    optparser.add_option('-f',"--force",
                         dest="force",
                         action="store_true",
                         help="force to run ClaRNA")

    optparser.add_option('-o',"--out_fn", type="string",
                         dest="out_fn",
                         default='inf.csv',
                         help="out csv file")

    (opts, args)=optparser.parse_args()

    if len(sys.argv) == 1:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)

    input_files = args[:] # opts.input_dir
    target_fn = opts.target_fn
    out_fn = opts.out_fn
    target_cl_fn = clarna_run(target_fn, opts.force)    

    f = open(out_fn, 'w')
    #t = 'target:' + os.path.basename(target_fn) + ' , rmsd_all\n'
    t = 'target,fn,inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC\n'
    print 'target, fn, inf_all, inf_stack, inf_WC, inf_nWC, SNS_WC, PPV_WC, SNS_nWC, PPV_nWC'
    f.write(t)
    for i in input_files:
        i_cl_fn = clarna_run(i, opts.force)
        scores = clarna_compare(target_cl_fn,i_cl_fn)
        print scores
        f.write(re.sub('\s+', ',', scores) + '\n')
    print 'csv was created! ', out_fn
