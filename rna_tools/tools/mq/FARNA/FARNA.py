#*-* coding: utf-8 *-*

"""Wrapper for ROSETTA software for structure prediction of small 
RNA sequences"""

import os, shutil, re
import subprocess

import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper, WrapperError
from rna_tools.tools.pdb_formatix.RosettaUtils import RosettaPDBFile
from rna_tools.tools.pdb_formatix.rebuilder import check_and_rebuild
from rna_tools.rna_tools_config import FARNA_PATH, FARNA_DB_PATH, FARNA_LORES


class FARNA(ProgramWrapper):
    """
    Wrapper class for running ROSETTA scoring function automatically.
    """
    program_name = 'farna'  
    src_bin = FARNA_PATH
    db_path = FARNA_DB_PATH
    input_fn = 'seq.fasta'
    input_file = ''
    best_energy = ''
    executable = 'rna_minimize'

    def __init__(self, sequence='test', seq_name='test', job_id=None):
        try:
            self.start_dir = os.getcwd()
        except OSError:        # directory was deleted or something like that
            pass
        super(FARNA, self).__init__(sequence, seq_name, job_id=job_id)

    def _prepare_stderr_stdout(self):
        # create output file
        self.output_file = os.path.join(self.path, 'stdout.txt')
        self.stdout = open(self.output_file, 'w')
        # create error file
        self.error_file = os.path.join(self.path, 'stderr.txt')
        self.stderr = open(self.error_file, 'w')

    def _prepare_files(self):
        # create input file
        self.input_file = open(os.path.join \
            (self.sandbox_dir, self.input_fn), 'w')\
                .write('>seq.fasta\n'+str(self.sequence).lower())
        self._prepare_stderr_stdout()
            
    def sandbox(self):
        #shutil.copytree(self.src_bin + os.sep + 'rosetta_source',
        #       self.sandbox_dir + os.sep + 'rosetta_source',
        #        symlinks=True)
        os.symlink(self.src_bin, self.sandbox_dir + os.sep + self.executable)
        #symlinks=True)

        #os.symlink(FARNA_DB_PATH,
        #        self.sandbox_dir + os.sep + 'rosetta_database')
        #os.system('chmod +x %s' % \
        #          os.path.join(self.sandbox_dir, self.executable))   
                               
    def run(self, pdb_file, hires, verbose=False, system=False):#, global_energy_score=True):
        """Compute FARNA potential for a single file

        Arguments:
            * pdb_file = path to pdb file
            * global_energy_score = True/False (See Output), default=True

        Output:
            * A list of energies, e.g::

              ['-21.721', '-0.899', '-20.961', '-84.498', '-16.574', '-180.939', '11.549', '7.475', '-17.257', '-306.324', '0.0', '0.0', '17.503', '0.0']

             ??? or a dictionary of lists of local scores, eg::

                {
                'N_BS': [17.0, -0.70039, -0.720981, -0.685238, -0.734146, ... ],
                'atom_pair_constraint': [0.0, -0.764688, -0.773833, ...],
                ...
                }

        """
        global_energy_score=True
        ftxt = open(pdb_file).read()
        ftxt = re.sub('TER\s+END\s+', 'TER', ftxt)
        ftxt = re.sub('END', 'TER', ftxt).strip()
        f = open(self.sandbox_dir + os.sep + 'tmp.pdb', 'w')
        f.write(ftxt)
        f.close()

        pdb_file = self.sandbox_dir + os.sep + 'tmp.pdb'

        if check_and_rebuild(pdb_file, self.sandbox_dir + os.sep + 'query.pdb'):
            self.pdb_fixes.append('rebuild_full_atom')
        pdb_file = RosettaPDBFile(pdb_path=self.sandbox_dir + os.sep + 'query.pdb')
        # get sequence from PDB file
        with open(self.sandbox_dir + os.sep + 'query.fasta', 'w') as f:
            f.write(pdb_file.get_fasta(lowercase=True))
        # create a ROSETTA ready PDB file
        pdb_file.make_rna_rosetta_ready()
        pdb_file.save(self.sandbox_dir + os.sep + 'query.pdb')
        self.pdb_fixes = pdb_file.fixes

        # run
        # os.chdir(self.sandbox_dir)
        self.flags = [self.sandbox_dir + os.sep + self.executable]

        # hires = True
        if hires == True: # False: # must be a string
            minimize_cmd = ' ' # -minimize_rna '
        else:
            minimize_cmd = ' -score:weights ' + FARNA_LORES + ' -minimize_rna '
            ## MM minimize_rna should be off or by option
            ## 2021 i'm not sure why? keep  -minimize_rna on here

        cmd = ' '.join([FARNA_PATH, '-constant_seed -database', self.db_path,
             minimize_cmd,
                        ' -ignore_zero_occupancy false ',
             '-s', self.sandbox_dir + os.sep + 'query.pdb',
             '-out:file:silent', self.sandbox_dir + os.sep + 'SCORE.out'])

        if verbose:
            print(cmd)
        self.log(cmd, 'debug')
        self.log('Running program')
        if system:
            os.system(cmd)
        else:
            out = subprocess.getoutput(cmd)

        self.log('Run finished')
        self.log(out, 'debug')

        self.get_result()

        #if global_energy_score: # ???
        results = []
        for i in list(self.result.keys()):
            results.append(str(self.result[i][0]))
        #return '\t'.join(results)
        return results
        
    def get_result(self):
        """Parse and get result from score file created during ROSETTA run

        All results are kept in self.result, but only global score is returned
        """
        f = open(self.sandbox_dir + os.sep + 'SCORE.out')
        output = f.read()
        f.close()
        lines = output.split('\n')
        lines = [l for l in lines if not l.startswith('REMARK')]
        # get names of different scores
        keys = lines[1].split()[1:-1]
        # get global scores
        global_scores = lines[2].split()[1:]
        # global scores are at index 0 in result, local are at 1--len(sequence)

        self.result = dict(list(zip(keys, [[float(s)] for s in global_scores[:len(keys)]])))
        
        ##for l in lines[3:-1]:
        #    scores_res = l.split()[2:-1]  # scores for a single residue
        #    for i in xrange(len(keys)):
        #        self.result[keys[i]].append(float(scores_res[i]))
        #return self.result['score'][0]
        return self.result

    def mqap(self, pdb):
        "Total weighted score:\s+(?P<ROSETTA_SCORE>[-\d.]+)"
        pass

    def cleanup(self):
        super(FARNA, self).cleanup()

# main
if __name__ == '__main__':
    fns = ['test.pdb']
    fns = ['1xjrA_M1.pdb', 'test.pdb'] 
    fns = ['3e5f_output4_01-000001_AA+ResnShift.pdb'] #2pcw_1_2chains.pdb']  # two chains
    for f in fns:
        f = 'test' + os.sep + f
        print('processing %s' % f)

        if 1:
            # mini false
            farna = FARNA()
            try:
                result = farna.run(f, False)
            except:
                result = 'error'
            print(result)
        if 1:
            # mini true
            farna = FARNA('', '')
            try:
                result = farna.run(f, True)
            except:
                result = 'error'
            print(result)

        #farna.cleanup()
