#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing eSCORE

Install::

   pip install barnaba


Output::

    $ baRNAba ESCORE --pdb ../test/1a9n.pdb --ff /Users/magnus/work/opt/barnaba/barnaba_201128/test/data/1S72.pdb
    # your output will be written to files with prefix outfile.ESCORE
    # KDE computed. Bandwidth=  0.25  using 10655 base-pairs# Loaded sample ../test/1a9n.pdb

    #     Frame       ESCORE
              0   4.1693e-01

The Escore could be also accessed via Python::

     from barnaba import escore
     Escore = escore.Escore([path_to_pdb])
     (..)
     # see example_12_escore.ipynb of barnaba package https://github.com/srnas/barnaba

"""
import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import RNA3DCNN_PATH, PYTHON3_PATH
from rna_tools.rna_tools_config import baRNAba_PATH

class eSCORE(ProgramWrapper):
    """
    Wrapper class for eSCORE
    """
    def __init__(self):
        super(eSCORE, self).__init__()

    def run(self, path_to_pdb, verbose=True):
        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        old_pwd = os.getcwd()
        #print baRNAba.parse() ## someday!

        os.chdir(self.sandbox_dir)
        
        self.log('eSCORE::start for %s' % self.sandbox_dir + '/query.pdb')
        # baRNAba ESCORE --pdb ../test/1a9n.pdb --ff /Users/magnus/work/opt/barnaba/barnaba_201128/test/data/1S72.pdb
        cmd = 'baRNAba ESCORE -o log.txt --ff ' + baRNAba_PATH + '/test/data/1S72.pdb --pdb ' + self.sandbox_dir + '/query.pdb &> /dev/null'
        os.system(cmd)
        self.log('eSCORE::Run finished')

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
