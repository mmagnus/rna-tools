#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions for computing RASP potential
"""
import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import RASP_PATH, RASP_PATH_DIR

# directory where this script is
# files for some commands are also here
DIRECTORY = os.path.dirname(__file__)


class RASP(ProgramWrapper):
    """
    Wrapper class for running RASP automatically.
    """
    program_name = ['rasp_fd', 'rasp_profile_fd']
    executable = ['rasp_fd', 'rasp_profile_fd']

    def __init__(self, job_id=None):
        super(RASP, self).__init__('', '', job_id=job_id, path=RASP_PATH)

    def run(self, path_to_pdb, global_energy_score=True, potentials=['c3', 'bb', 'bbr', 'all'], handler=True, verbose=False):
        """Compute RASP potential for a single file

        ``rasp_fd`` generates: 
            ``output_all.txt``::

            -9869.61   66447       -0.148534       0       0       0

        ``rasp_profile_fd`` generates:

            ``profile_all.txt``::

                C	    1	    R	    -775.038
                C	    2	    R	    -1164.22
                U	    3  	    R	    -2054.17
                G           4	    R	    -1601.13
                [..]
        
        Input:
            * path_to_pdb = path to PDB file
            * potential_type = all, bbr or c3
            * global_energy_score = True/False (See Output), default=True
            * hander = True (if you use PYRO)/FALSE otherwise

        Output:
            * the output depends on global_enery_score value T/F
              You might get:

              ``-9869.61``

              ``profile [('C', 1, 'R', -775.03800000000001), ('C', 2, 'R', -1164.22) ..``

        """
        # step 0
        if handler:
            copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        else:  # if Pyro needs this because it cannot access file on client
            # TODO: make a function decorator with this
            f = open(self.sandbox_dir + os.sep + 'query.pdb', 'w')
            f.write(path_to_pdb)
            f.close()

        old_pwd = os.getcwd()
        os.chdir(self.sandbox_dir)
        pdb_file = PDBFile(pdb_path=os.path.join(self.sandbox_dir, 'query.pdb'))
        pdb_file.resname_check_and_3to1()
        pdb_file.save(os.path.join(self.sandbox_dir, 'query.pdb'))
        self.pdb_fixes = pdb_file.fixes

        scores = [] # colect all scores!

        if global_energy_score:
            # : 
            # this order is important, should not be changed!!! (based on this order mqaprnadb_load_from.csv reads the data
            for potential_type in potentials:
                self.log('Running "global_energy_score"')
                if run_command(RASP_PATH + os.sep + 'bin' + os.sep + self.executable[0],
                            ['-e', potential_type, '-p', 'query.pdb'],
                            env={'RASP': RASP_PATH},
                               stdout_file='output_' + potential_type + '.txt', verbose=verbose):
                    print('Computing global energy failed')
                    return -1
                else:
                    f = open('output_' + potential_type + '.txt')
                    result = [l.split() for l in f.readlines()]
                    scores += result[0] # ['-1.9657', '144', '-0.0136507', '0.794594', '1.10601', '-2.49573', '-2309.53', '8595', '-0.268707', '-2199.69', '18.6549', '-5.88841', '-5055.01', '19931', '-0.253625', '-4665.96', '44.262', '-8.7896', '-9869.61', '66447', '-0.148534', '0', '0', '0']
                    #result = float(result[0][2])
                    f.close()

        else:
            self.log('Running profile calculation')
            if run_command(self.sandbox_dir + os.sep + self.executable[1],
                           ['-e', potential_type, '-p', 'query.pdb'],
                           env={'RASP': RASP_PATH},
                           stdout_file='profile_' + potential_type + '.txt'):
                print('Computing profile failed')
                return -1
            else:
                f = open('profile_' + potential_type + '.txt')
                result = [l.split() for l in f.readlines()]
                try:
                    result = [(int(i[1]), float(i[3])) for i in result]
                except IndexError:
                    result = [(int(i[1]), float(i[2])) for i in result]
                f.close()

        self.log('Run finished')
        os.chdir(old_pwd)
        #return '\t'.join(scores)
        return scores
    #return result or -1  # return -1 if result is 0 or None

    def colour_by_local_score(self, path_to_pdb, potential_type='all'):
        local_scores = self.run(path_to_pdb, potential_type, global_energy_score=False)
        local_scores = [float(i[1]) for i in local_scores]
        with open(path_to_pdb) as f:
            pdb_file = PDBFile(pdb_string=f.read())
        pdb_file.set_residues_bfactor(local_scores)
        return pdb_file.pdb_string
