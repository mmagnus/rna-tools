#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
from shutil import copyfile, copytree, rmtree
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import rsRNASP_PATH

import logging

# directory where this script is
# files for some commands are also here
#try:

#except NameError:
DIRECTORY = os.path.dirname(__file__)

class rsRNASP(ProgramWrapper):
    """
    Wrapper class for running rsRNASP
    """
    program_name = 'rsRNASP'

    def __init__(self, job_id=None):
        super(rsRNASP, self).__init__('', '', job_id=job_id)

    def run(self, path_to_pdb, verbose=1):
        self.logger.setLevel(logging.INFO)

        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        copyfile(rsRNASP_PATH + '/rsRNASP', self.sandbox_dir + os.sep + 'rsRNASP')
        os.symlink(rsRNASP_PATH + '/data', self.sandbox_dir + os.sep + 'data')
        os.chmod(self.sandbox_dir + os.sep + 'rsRNASP', 0o777)
        old_pwd = os.getcwd()

        #pdb_file = PDBFile(pdb_path=os.path.join(self.sandbox_dir, 'query.pdb'))
        #pdb_file.resname_check_and_3to1()
        #pdb_file.save(os.path.join(self.sandbox_dir, 'query.pdb'))
        #self.pdb_fixes = pdb_file.fixes
        os.chdir(self.sandbox_dir)

        err = ''
        out = ''
        self.log('rsRNASP::start for %s' % self.sandbox_dir + '/query.pdb')
        os.chdir(self.sandbox_dir)
        #cmd = 'export rsRNASP=' + rsRNASP_PATH + ';'
        # hack for m1 running arch -x86_64 rsRNASP
        cmd = 'cd ' + rsRNASP_PATH + ' && ./rsRNASP ' + self.sandbox_dir + '/query.pdb' + ' ' + self.sandbox_dir + '/output.txt >> ' + self.sandbox_dir + '/std.txt'
        os.system(cmd)
        #    stdout_file=self.sandbox_dir + '/std.txt', stderr_file=self.sandbox_dir + '/err.txt', verbose=verbose):
        #    err = open(self.sandbox_dir + '/err.txt').read().strip()
        try:
            f = open(self.sandbox_dir + '/output.txt')
        except:
            return -1
        try:
            result = float(f.read().split()[1])  # 'query.pdb     -3055.902390\n'
        except:
            result = -1
        if err:
            self.log('rsRNASP::err %s' % err)
        self.log('rsRNASP::result %s' % str(result))
        self.log('rsRNASP::Run finished')

        rmtree(self.sandbox_dir)
        os.chdir(old_pwd)
        return result # return -1 if result is 0 or None

def main():
    wrapper = rsRNASP()
    try:
        result = wrapper.run('../test' + os.sep + '1a9n.pdb')
        if result:
            print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

if '__main__' == __name__:
    main()
