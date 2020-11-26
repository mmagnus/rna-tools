#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing Dfire potential
"""
import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import dfire_PATH

class Dfire(ProgramWrapper):
    """
    Wrapper class for Dfire.
    """
    max_seq_len = 100000  # I don't know about any restriction

    def __init__(self, sequence, seq_name, job_id=None):
        super(Dfire, self).__init__(sequence, seq_name, job_id=job_id)

    def run(self, path_to_pdb, verbose=False):
        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        old_pwd = os.getcwd()
        os.chdir(self.sandbox_dir)
        # self.log('dfire::start for %s' % self.sandbox_dir + '/query.pdb')
        cmd = "export DFIRE_RNA_HOME=" + dfire_PATH + "; "
        cmd += dfire_PATH + '/bin/DFIRE_RNA ' + self.sandbox_dir + '/query.pdb >'  + self.sandbox_dir + '/log.txt.dfire.rna'
        if verbose: print(cmd)
        os.system(cmd)

        # self.log('dfire::Run finished')

        for line in open(self.sandbox_dir + '/log.txt.dfire.rna'):
            score = line.strip().split()[1] # /tmp/tmplDNztf/query.pdb -12480.188898
        if verbose: print(open(self.sandbox_dir + '/log.txt.dfire.rna').read())
        os.chdir(old_pwd)
        return float(score)


def main():
    wrapper = Dfire('', '')
    try:
        result = wrapper.run('test' + os.sep + '1a9n.pdb')
        if result:
            print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

    print((wrapper.run('test' + os.sep + '1a9nR.pdb')))
    print((wrapper.run('test' + os.sep + '3b58ABC.pdb')))
    print((wrapper.run('test' + os.sep + '5e3hBC.pdb')))
    print((wrapper.run('test' + os.sep + 'S_000001_000.pdb')))
    

if '__main__' == __name__:
    main()
