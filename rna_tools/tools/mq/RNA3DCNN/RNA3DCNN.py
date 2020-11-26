#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing RASP potential
"""
import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import RNA3DCNN_PATH, PYTHON3_PATH


class RNA3DCNN(ProgramWrapper):
    """
    Wrapper class for RNA3DCNN.
    """
    max_seq_len = 100000  # I don't know about any restriction

    def __init__(self, sequence, seq_name, job_id=None):
        super(RNA3DCNN, self).__init__(sequence, seq_name, job_id=job_id)

    def run(self, path_to_pdb, verbose=False):
        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        old_pwd = os.getcwd()

        os.chdir(self.sandbox_dir)

        self.log('start for %s' % self.sandbox_dir + '/query.pdb', level="debug")

        # 4. To print scores of each nucleotide and total scores, use flag "-local 1"
        # 5. To print only total scores, use flag "-local 0"
        # For example:<br />
        # python Main.py -pl pdblist -model RNA3DCNN_MD.hdf5 -local 0<br />
        cmd = PYTHON3_PATH + ' ' + RNA3DCNN_PATH + '/Main.py ' + \
        ' -pn ' + self.sandbox_dir + '/query.pdb ' + \
        ' -model ' + RNA3DCNN_PATH + '/RNA3DCNN_MD.hdf5 ' + \
        ' -local 0 2>>  '  + self.sandbox_dir + '/log.txt >>' + self.sandbox_dir + '/log.txt'
        if verbose:
            print(cmd)
        os.system(cmd)

        """
        Total params: 4,282,801
        Trainable params: 4,282,801
        Non-trainable params: 0
        _________________________________________________________________
        Total score for /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpO_jRVR/query.pdb is  6.3262305
        """
        self.log('Run finished')
        score = open(self.sandbox_dir + '/log.txt').read().split()[-1]
        if verbose:
            print(score)
        os.chdir(old_pwd)
        try:
            return float(score)
        except:
            error = 'Error: problem with the file'
            with open(self.sandbox_dir + '/log.txt') as f:
                error += f.read()
            return error

def main():
    wrapper = RNA3DCNN('', '')
    try:
        result = wrapper.run('test' + os.sep + '1a9n.pdb')
        if result:
            print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

if '__main__' == __name__:
    main()
