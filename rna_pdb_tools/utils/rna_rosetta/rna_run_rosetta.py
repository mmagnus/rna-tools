"""**run_rosetta** - wrapper to ROSETTA tools for RNA modeling

Based on C. Y. Cheng, F. C. Chou, and R. Das, Modeling complex RNA tertiary folds with Rosetta, 1st ed., vol. 553. Elsevier Inc., 2015.
http://www.sciencedirect.com/science/article/pii/S0076687914000524 """

import argparse
import os
import glob
import subprocess

ROOT_DIR_MODELING = '/home/magnus/'

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='file')
    return parser

def prepare_folder(header,seq,ss,path):
    """Make folder for you job, with seq.fa, ss.fa"""
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

def make_helices(path):
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

def run_rosetta(path):
    """
    change to 50 per job (!) change manually right now
    50 * 100 = 10k ?
    """
    # get list line
    helices = open('CMDLINES').readlines()[-1].replace('#','')
    cmd = 'rna_denovo_setup.py -fasta seq.fa -secstruct_file ss.fa -cycles 20000 ' + helices
    print cmd
    os.system(cmd)
    # change to 50 per job (!)
    # 50 * 100 = 10k ?
    cmd = 'rosetta_submit.py README_FARFAR o 200 72 ' + header[:3]
    print cmd
    os.system(cmd)    

    os.system('chmod +x ./qsubMINI')
    os.system('./qsubMINI &')

def run():
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
        os.chdir(path)
        
        #prepare_folder(header, seq, ss, path)
        #make_helices(path)
        run_rosetta(path)

        os.chdir(curr)

#main
if __name__ == '__main__':
    run()
