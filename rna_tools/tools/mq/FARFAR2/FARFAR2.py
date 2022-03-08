#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper for ROSETTA/FARFAR2 (2020 update) software for structure prediction
of small RNA sequences.

Scores::

    ------------------------------------------------------------
     Scores                       Weight   Raw Score Wghtd.Score
    ------------------------------------------------------------
     fa_atr                       0.210   -1081.212    -227.054
     fa_rep                       0.200     180.360      36.072
     fa_intra_rep                 0.003    1504.776       4.364
     lk_nonpolar                  0.250     -26.712      -6.678
     fa_elec_rna_phos_phos        1.700      -3.955      -6.724
     rna_torsion                  1.000      63.080      63.080
     suiteness_bonus              1.000       0.000       0.000
     rna_sugar_close              0.820      13.043      10.696
     fa_stack                     0.130   -1144.488    -148.783
     stack_elec                   0.760     -17.443     -13.257
     geom_sol_fast                0.170     355.629      60.457
     hbond_sr_bb_sc               0.960     -12.355     -11.861
     hbond_lr_bb_sc               0.960      -0.562      -0.539
     hbond_sc                     0.960     -76.482     -73.423
     ref                          1.000     257.400     257.400
     free_suite                   2.000       0.000       0.000
     free_2HOprime                1.000       0.000       0.000
     intermol                     1.000       0.000       0.000
     other_pose                   1.000       0.000       0.000
     loop_close                   1.000       0.000       0.000
     linear_chainbreak            5.000       0.000       0.000
    ---------------------------------------------------
     Total weighted score:                      -56.252
    protocols.rna.denovo.movers.RNA_Minimizer: RNA minimizer finished in 230 seconds.

csv::

    rna_csv.py --flat FARFAR2_hires-3.7.19+5.gd7c6fd7.dirty-mx.local.csv
    id 1
    fn 1y26X_output4_01-000085_AA.pdb
    ff2_score_hires -56.252
    ff2_fa_atr -227.054
    ff2_fa_rep 36.072
    ff2_fa_intra_rep 4.364
    ff2_lk_nonpolar -6.678
    ff2_fa_elec_rna_phos_phos -6.724
    ff2_rna_torsion 63.08
    ff2_suiteness_bonus 0.0
    ff2_rna_sugar_close 10.696
    ff2_fa_stack -148.783
    ff2_stack_elec -13.257
    ff2_geom_sol_fast 60.457
    ff2_bond_sr_bb_sc -11.861
    ff2_hbond_lr_bb_sc -0.539
    ff2_hbond_sc -73.423
    ff2_ref 257.4
    ff2_free_suite 0.0
    ff2_free_2HOprime 0.0
    ff2_intermol 0.0
    ff2_other_pose 0.0
    ff2_loop_close 0.0
    ff2_linear_chainbreak_hires 0.0

"""

import os, shutil, re
import subprocess

import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper, WrapperError
from rna_tools.tools.pdb_formatix.RosettaUtils import RosettaPDBFile
from rna_tools.tools.pdb_formatix.rebuilder import check_and_rebuild
from rna_tools.rna_tools_config import FARFAR2_PATH, FARFAR2_DB_PATH, FARFAR2_LORES


class FARFAR2(ProgramWrapper): 
    """
    Wrapper class for running ROSETTA scoring function automatically.
    """
    program_name = 'farna'  
    src_bin = FARFAR2_PATH
    db_path = FARFAR2_DB_PATH
    input_fn = 'seq.fasta'
    input_file = ''
    best_energy = ''
    executable = FARFAR2_PATH

    def __init__(self, job_id=None):
        try:
            self.start_dir = os.getcwd()
        except OSError:        # directory was deleted or something like that
            pass
        super(FARFAR2, self).__init__('', '', job_id=job_id)

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
        #os.symlink(self.src_bin, self.sandbox_dir + os.sep + self.executable)
        #symlinks=True)

        #os.symlink(FARNA_DB_PATH,
        #        self.sandbox_dir + os.sep + 'rosetta_database')
        #os.system('chmod +x %s' % \
        #          os.path.join(self.sandbox_dir, self.executable))
        pass
                               
    def run(self, pdb_file, hires, verbose=False):#, global_energy_score=True):
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
        #os.chdir(self.sandbox_dir)
        self.flags = [self.sandbox_dir + os.sep + self.executable]

        if hires == False: # must be a string
            minimize_cmd = ' -score:weights ' + FARFAR2_LORES + ' ' # # -minimize_rna' ##
            # MM minimize_rna should be off or by option
        else:
            minimize_cmd = ' '
       
        # self.sandbox_dir + os.sep +
        # verbose = True
        cmd = ' '.join([self.executable , '-database', self.db_path,
             minimize_cmd,
                        ' -ignore_zero_occupancy false -constant_seed ',
             '-s', self.sandbox_dir + os.sep + 'query.pdb',
             '-out:file:silent', self.sandbox_dir + os.sep + 'SCORE.out'])

        self.log(cmd, 'debug')
        if verbose: print(cmd)
        self.log('Running program')

        system = False # for debuggings
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
            results.append(float(self.result[i][0]))  # float vs string
        #return '\t'.join(results)
        return results
        
    def get_result(self):
        """Parse and get result from score file created during ROSETTA run

        All results are kept in self.result, but only global score is returned

        """
        try:
            with open(self.sandbox_dir + os.sep + 'SCORE.out') as f:
                lines = f.readlines() #().split('\n')
        except:
            raise Exception('Problem with gettting the results, run with --verbose')
        lines = [l for l in lines if not l.startswith('REMARK')]
        #print(lines)
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
        super(FARFAR2, self).cleanup()

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
            farna = FARFAR2()
            result = farna.run(f, False, verbose=True)
            print(result)

        if 1:
            # mini true
            farna = FARFAR2()
            result = farna.run(f, True, verbose=True)
            print(result)

        #farna.cleanup()
