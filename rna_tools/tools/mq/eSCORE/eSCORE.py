#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing eSCORE

Install::

   pip install barnaba # tested with barnaba==0.1.7

Output::

    $ baRNAba ESCORE --pdb ../test/1a9n.pdb --ff /Users/magnus/work/opt/barnaba/barnaba_201128/test/data/1S72.pdb
    # your output will be written to files with prefix outfile.ESCORE
    # KDE computed. Bandwidth=  0.25  using 10655 base-pairs# Loaded sample ../test/1a9n.pdb

    #     Frame       ESCORE
              0   4.1693e-01

The eSCORE could be also accessed via Python::

     from barnaba import escore
     Escore = escore.Escore([path_to_pdb])
     (..)
     # see example_12_escore.ipynb of barnaba package https://github.com/srnas/barnaba

It does not work on M1 mac (problem to compile mdtraj).
"""

import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import baRNAba_data_PATH

import subprocess

def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


class eSCORE(ProgramWrapper):
    """
    Wrapper class for eSCORE
    """
    def __init__(self):
        super(eSCORE, self).__init__()

    def run(self, path_to_pdb, verbose=False):
        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        old_pwd = os.getcwd()
        #print baRNAba.parse() ## someday!

        os.chdir(self.sandbox_dir)

        self.log('eSCORE::start for %s' % self.sandbox_dir + '/query.pdb', verbose=verbose)
        # baRNAba ESCORE --pdb ../test/1a9n.pdb --ff /Users/magnus/work/opt/barnaba/barnaba_201128/test/data/1S72.pdb
        #TOFIX!
        baRNAba_data_PATH='/home/mqapRNA/mqaprna_env/opt/barnaba/examples/DATA/1S72.pdb'
        cmd = 'barnaba ESCORE -o log.txt --ff ' +  baRNAba_data_PATH + ' --pdb ' + self.sandbox_dir + '/query.pdb &> /dev/null'
        if verbose:
            print(cmd)
        exe(cmd)

        self.log(cmd, verbose=verbose)
        self.log('eSCORE::Run finished', verbose=verbose)

        for line in open(self.sandbox_dir + '/log.txt.ESCORE.out'):
            if not line.startswith('#'):
                score = line.strip().split()[1] # 0   4.1693e-01 # 0.41693
        os.chdir(old_pwd)
        return float(score)

def main():
    wrapper = eSCORE()
    try:
        result = wrapper.run(os.path.abspath('../test' + os.sep + '1a9n.pdb'))
        print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

if '__main__' == __name__:
    main()
