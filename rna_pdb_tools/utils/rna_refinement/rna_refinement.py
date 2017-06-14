#!/usr/bin/python
# -*- coding: utf-8 -*-

"""RNA refinement with QRNAS

Install (http://genesilico.pl/qrnas)

Be default the script searches QRNAS in (rna-pdb-tools)/opt/qrnas/ .

Usage of QRNA::

    QRNA - Quick Refinement of Nucleic Acids (0.2 alpha)
         by Juliusz Stasiewicz (jstasiewicz@genesilico.pl)

    To use type:
      QRNA -i <input PDBfile> [-o <output PDBfile>] [-c <configfile>] [-p] [-m <restraintsfile>]
    OR specify <input PDBfile>, <output PDBfile> and <restraintsfile> in <configfile> and type just:
      QRNA -c <configfile>

TODO:

- clean up the output structure
- configuration should not be hardcoded

"""

from __future__ import print_function
import argparse
import re
import os
import subprocess
import random
import string

from shutil import copyfile

PATH = "/Users/magnus/work/src/rna-pdb-tools/"
QRNAS_PATH = PATH + '/opt/qrnas/'

class QRNAS:
    """QRNAS"""
    def run(self, inputfile, outputfile, steps = 10):
        """pdb_txt - content of pdb file,
        steps - nr of steps of simulation, 500 by default.

        NSTEPS must be uncommented."""
        cwd = os.getcwd()
        # get config
        conftxt = open(QRNAS_PATH + os.sep + 'configfile.txt').read()
        conftxt_tmp = re.sub('NSTEPS.+\d+', 'NSTEPS   ' + str(steps), conftxt)

        JOB_ID = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(20))
        JOB_PATH = QRNAS_PATH + os.sep + 'jobs/' + JOB_ID
        os.makedirs(JOB_PATH)

        # get temp config
        with open(JOB_PATH + os.sep + 'configfile.txt','w') as f:
            f.write(conftxt_tmp)

        # copy input to qrnas folder
        qrnas_inputfile = QRNAS_PATH + os.sep + 'jobs/' + JOB_ID + os.sep + os.path.basename(inputfile)
        copyfile(inputfile, qrnas_inputfile)

        os.chdir(QRNAS_PATH)
        cmd = './QRNAS -i jobs/' + JOB_ID + os.sep + os.path.basename(inputfile) + \
          ' -c jobs/' + JOB_ID + os.sep + 'configfile.txt ' + \
          ' -o jobs/' + JOB_ID + os.sep + os.path.basename(inputfile).replace('.pdb', '.refi.pdb')
        print(cmd)
        subprocess.call(cmd, shell=True)

        os.chdir(cwd)
        copyfile(QRNAS_PATH + '/jobs/' + JOB_ID + os.sep + os.path.basename(inputfile).replace('.pdb', '.refi.pdb'), outputfile)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', help="input pdb file")     
    parser.add_argument('outputfile', help="output pdb file")
    return parser
#main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    q = QRNAS()
    q.run(args.inputfile, args.outputfile)
