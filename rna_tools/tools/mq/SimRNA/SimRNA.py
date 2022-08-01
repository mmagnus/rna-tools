#!usr/bin python
#-*- coding: utf-8 -*-
"""Module for getting SimRNA energy for an RNA molecule.
"""
import os
import types
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import SIMRNA_PATH, SIMRNA_DATA_PATH
from rna_tools.rna_tools_lib import RNAStructure

class SimRNA(ProgramWrapper):
    """
    Wrapper class for running SimRNA automatically.
    """
    program_name = 'SimRNA'
    max_seq_len = 2000
    executable = None
    src_bin = SIMRNA_PATH + os.sep + program_name
    frame = ''

    def __init__(self, job_id=None):
        super(SimRNA, self).__init__('', '', job_id=job_id)

    def _prepare_files(self):
        try:
            self.log(SIMRNA_DATA_PATH + '->' + self.sandbox_dir + os.sep + 'data')
            os.symlink(SIMRNA_DATA_PATH, self.sandbox_dir + os.sep + 'data')
        # go on if symlink is already there
        except OSError:
            self.log('Symlink already exists', 'debug')

    def run(self, pdb_file, numSteps, verbose=False):
        """Compute SimRNA energy for a single file

        Arguments:
            * pdb_file = path to a PDB file
            * numSteps, number of SimRNA steps
            
        Output:
            * global energy as a floating point number
        """
        self.log('Running program')
        self.log('input - pdb_file: %s' % pdb_file)

        #copyfile(pdb_file, self.sandbox_dir + os.sep + 'query.pdb')
        s = RNAStructure(pdb_file)
        s.get_rnapuzzle_ready()
        s.write(self.sandbox_dir + os.sep + 'query.pdb', False)

        try:
            old_pwd = os.getcwd()
        except OSError:
            old_pwd = None
        os.chdir(self.sandbox_dir)

        if verbose: print('simrna::', self.sandbox_dir)

        # This is the place why you run SimRNA (self.src_bin)
        # you can add options like '-c', '~/tmp/config' 
        # following flags mean:
        # -n 0 = 0 steps of simulation
        # -R 0 = set random seed to a constant,
        # results should be the same every time
        if run_command(self.src_bin,
                    ['-p', self.sandbox_dir + os.sep + 'query.pdb',
                        '-n', str(numSteps), '-R', '0'], # 16000
                    stderr_file=self.sandbox_dir + os.sep + 'output.txt'):
            self.log('Run failed', 'error')
            
            # the place where I read the output file (stderr_file) of SimRNA
            f = open(self.sandbox_dir + os.sep + 'output.txt')
            stderr = f.read()
            f.close()
            self.log('STDERR:\n' + stderr, 'debug')
            if old_pwd:
                os.chdir(old_pwd)
            return -1
        if old_pwd:
            os.chdir(old_pwd)
        self.log('Run finished')
        
        # parse output
        scores = []
        f = open(self.sandbox_dir + os.sep + 'output.txt')
        lines_temp = f.read().split('\n')
        lines = lines_temp
        
        if verbose: print('simrna::', '\n'.join(lines))

        # the output is dicombosed into score
        total_enrg = [l for l in lines if 'Total energy of system' in l][-1].split(' ')[-1]
        bb = [l for l in lines if 'Base-Base interactions energy' in l][-1].split(' ')[-1]
        ss = [l for l in lines if 'where:  short stacking energy' in l][-1].strip().split(' ')[-1]
        bbc = [l for l in lines if 'Base-Backbone interact. energy' in l][-1].strip().split(' ')[-1]
        lge = [l for l in lines if 'Local geometry energy' in l][-1].strip().split(' ')[-1]

        # for VM
        if 1: # diff between SimRNA versions?
            bcp = [l for l in lines if "bonds (distance)    C4'-P" in l][-1].strip().split(' ')[-1]
            bpc = [l for l in lines if "bonds (distance)    P-C4'" in l][-1].strip().split(' ')[-1]
            fa = [l for l in lines if "flat angles     C4'-P-C4'" in l][-1].strip().split(' ')[-1]
            fa2 = [l for l in lines if "flat angles      P-C4'-P" in l][-1].strip().split(' ')[-1]
            tor = [l for l in lines if "torsion eta vs torsion theta" in l][-1].strip().split(' ')[-1]
        # for MAC
        else:
            bcp = [l for l in lines if "bonds (distance)  C4'-P" in l][-1].strip().split(' ')[-1]
            bpc = [l for l in lines if "bonds (distance)  P-C4'" in l][-1].strip().split(' ')[-1]
            fa = [l for l in lines if "flat angles   C4'-P-C4'" in l][-1].strip().split(' ')[-1]
            fa2 = [l for l in lines if "flat angles    P-C4'-P" in l][-1].strip().split(' ')[-1]
            tor = [l for l in lines if "tors. eta vs tors. theta" in l][-1].strip().split(' ')[-1]

        lim = '' # empty at the moment '' # [l for l in lines if "limit. sphere exceed penalty:" in l][-1].strip().split(' ')[-1]
        cc =[l for l in lines if "current chain energy:" in l][-1].strip().split(' ')[-1]
        scores += [numSteps] + [total_enrg] + [bb] + [ss] + [bbc] + [lge] + [bcp] + [bpc] + [fa] + [fa2] + [tor] + [lim] + [cc]
        f.close()
        if verbose: print(' simrna::steps:', numSteps)
        #s = '\t'.join(scores)
        #if verbose: print ' simrna::scores:', s
        #return s
        return scores

    def cleanup(self):
        super(SimRNA, self).cleanup()


if __name__ == '__main__':
    wrapper = SimRNA()
    result = wrapper.run('../test/1xjrA_M1.pdb', 0, verbose=False) # or aa
    print(result)
