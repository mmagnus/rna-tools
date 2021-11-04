#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing 3dRNAscore

Install::

    make # make clean if you don't have clean

At Mac there is a problem (201128, 211103)::

    (py37) [mx] example$ ../bin/3dRNAscore -s:l score.in >score.txt
    [1]    69818 segmentation fault  ../bin/3dRNAscore -s:l score.in > score.txt

"""
import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import RNAscore_PATH

# directory where this script is
# files for some commands are also here
#try:

#except NameError:
DIRECTORY = os.path.dirname(__file__)

class RNAscore(ProgramWrapper):
    """
    Wrapper class for running 3dRNAscore
    """
    program_name = '/bin/3dRNAscore'

    def __init__(self, job_id=None):
        super(RNAscore, self).__init__('', '', job_id=job_id)

    def run(self, path_to_pdb, verbose=1):
        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        old_pwd = os.getcwd()

        #pdb_file = PDBFile(pdb_path=os.path.join(self.sandbox_dir, 'query.pdb'))
        #pdb_file.resname_check_and_3to1()
        #pdb_file.save(os.path.join(self.sandbox_dir, 'query.pdb'))
        #self.pdb_fixes = pdb_file.fixes
        os.chdir(self.sandbox_dir)

        err = ''
        out = ''
        self.log('3dRNAscore::start for %s' % self.sandbox_dir + '/query.pdb')
        if run_command(RNAscore_PATH + os.sep + self.program_name,
                        ['-s', self.sandbox_dir + '/query.pdb'],
                        env={'RNAscore': RNAscore_PATH},
                       stdout_file=self.sandbox_dir + '/output.txt', stderr_file=self.sandbox_dir + '/err.txt', verbose=verbose):
            err = open(self.sandbox_dir + '/err.txt').read().strip()
            print('3dRNAscore::err', err)
            result = -1
        else:
                f = open(self.sandbox_dir + '/output.txt')
                result = float(f.read())
                f.close()

        self.log('3dRNAscore::err %s' % err)
        self.log('3dRNAscore::result %s' % str(result))
        self.log('3dRNAscore::Run finished')

        os.chdir(old_pwd)
        return result # return -1 if result is 0 or None

def main():
    wrapper = RNAscore()
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
