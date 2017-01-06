#!/usr/bin/env python

"""**run_rosetta** - wrapper to ROSETTA tools for RNA modeling

Based on C. Y. Cheng, F. C. Chou, and R. Das, Modeling complex RNA tertiary folds with Rosetta, 1st ed., vol. 553. Elsevier Inc., 2015.
http://www.sciencedirect.com/science/article/pii/S0076687914000524

The script makes (1) a folder for you job, with seq.fa, ss.fa, input file is copied as input.fa to the folder (2) make helices (3) prepare rosetta input files (4) sends jobs to the cluster."""

import argparse
import textwrap
import os
import glob
import subprocess
import shutil

ROOT_DIR_MODELING = '/home/magnus/'

def get_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--init', help='prepare _folder with seq and ss',
                        action='store_true')
    parser.add_argument('-e', '--helices', help='produce h(E)lices',
                        action='store_true')
    parser.add_argument('-r', '--rosetta', help='prepare rosetta files (still you need `go` to send jobs to a cluster)',
                        action='store_true')
    parser.add_argument('-g', '--go', help= go.__doc__,
                        action='store_true')

    parser.add_argument('file', help= textwrap.dedent("""file:\n>a04\nUAUAACAUAUAAUUUUGACAAUAUGGGUCAUAAGUUUCUACCGGAAUACCGUAAAUAUUCUGACUAUGUAUA\n((((.((((...((.((((.......)))).))........(.(((((.......))))).)..))))))))"""))
    return parser

def prepare_folder(args,header,seq,ss,path):
    """Make folder for you job, with seq.fa, ss.fa, input file is copied as input.fa to the folder"""
    d = path
    try:
        os.mkdir(d)
        print d, 'created'
    except OSError:
        print d, 'created is already created'
        pass

    with open(d + "seq.fa","w") as f:
        f.write(header +'\n')
        f.write(seq)

    with open(d + "ss.fa","w") as f:
        f.write(ss)
    print 'Seq & ss created'
    shutil.copyfile(args.file, d + 'input.fa')

def prepare_helices():
    """Make helices (wrapper around 'helix_preassemble_setup.py')"""
    cmd = 'helix_preassemble_setup.py -secstruct ss.fa -fasta seq.fa'
    os.system(cmd)
    #
    helix_runs = glob.glob('*RUN')
    print helix_runs

    f = open('HRUNS', 'w')
    for h in helix_runs:
        f.write(open(h).read().strip() + ' & \n')
    f.close()
    
    p = subprocess.Popen(open('HRUNS').read(), shell=True, stderr=subprocess.PIPE)
    p.wait()
    stderr = p.stderr.read().strip()
    if stderr:
        print stderr

def prepare_rosetta(header):
    """Repare ROSETTA using rna_denovo_setup.py

    change to 50 per job (!) change manually right now
    50 * 100 = 10k ?
    40 * 500 = 20k
    """
    # get list line
    helices = open('CMDLINES').readlines()[-1].replace('#','')
    nstruct = 40 #500
    njobs = 500
    cmd = 'rna_denovo_setup.py -fasta seq.fa -secstruct_file ss.fa -cycles 20000 -no_minimize -nstruct ' + str(nstruct) + ' ' + helices
    print cmd
    os.system(cmd)
    # change to 50 per job (!)
    # 50 * 100 = 10k ?
    cmd = 'rosetta_submit.py README_FARFAR o ' + str(njobs) + ' 100 ' + header[:3]
    print cmd
    os.system(cmd)    

def go():
    """send jobs to a cluster (run qsubMINI)"""
    os.system('chmod +x ./qsubMINI')
    os.system('./qsubMINI')
    
def main():
    """Pipline for modeling RNA"""
    args = get_parser().parse_args()
    
    if args.file:
        f = open(args.file)
        header = f.readline().strip().replace('>','')
        seq = f.readline().strip()
        ss = f.readline().strip()

        print 'run rosetta for:'
        print header
        print seq
        print ss

        path = ROOT_DIR_MODELING + os.sep + header + os.sep
        curr = os.getcwd()
        if args.init:
            prepare_folder(args,header, seq, ss, path)
        try:
            os.chdir(path)
        except OSError:
            print 'You have to make a folder first! use --init'

        if args.helices:
            prepare_helices()
        if args.rosetta:
            prepare_rosetta(header)
        if args.go:
            go()

        os.chdir(curr)
#main
if __name__ == '__main__':
    main()
