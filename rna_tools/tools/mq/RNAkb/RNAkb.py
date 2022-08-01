#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing RNAkb potential

It seems that this is impossible to run RNAkb in full atom mode. So this works only in 5 pt (5 points/atom per residue) mode.

https://gromacs.bioexcel.eu/t/fatal-error-an-input-file-contains-a-line-longer-than-4095-characters/1397/2
"""
import os
import sys
import re

from math import isnan
from shutil import copyfile

from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.tools.md.rnakb_utils import make_rna_rnakb_ready, fix_gromacs_gro, prepare_groups, format_score_mdp
from rna_tools.rna_tools_config import GRMLIB_PATH, DB_AA_PATH, DB_5PT_PATH, WRAPPERS_PATH, GROMACS_LD_PATH, GROMACS_PATH_BIN, TMP_PATH

import tempfile

class RNAkb(ProgramWrapper):
    """Wrapper class for running RNAkb automatically.
    """
    executable = ['/usr/local/gromacs/bin/' + exe for exe in ['pdb2gmx', 'make_ndx', 'editconf', 'grompp', 'mdrun']]

    def __init__(self, job_id='',  sandbox=False):
        SANDBOX_PATH = tempfile.TemporaryDirectory().name # '.' or './sandbox' for sandbox here
        super(RNAkb, self).__init__('', '', job_id=job_id)

        # sandbox q-&-d
        if sandbox:
            try:
                os.mkdir(SANDBOX_PATH)
            except:
                pass
            # os.system('rm ' + SANDBOX_PATH + '/*')
            self.sandbox_dir= SANDBOX_PATH

    def log_stdout_stderr(self):
        stdout_file = open(self.sandbox_dir + os.sep + 'out.txt')
        stderr_file = open(self.sandbox_dir + os.sep + 'err.txt')
        stdout = stdout_file.read()
        stderr = stderr_file.read()
        self.log('STDOUT:\n' + stdout, 'debug')
        self.log('STDERR:\n' + stderr, 'debug')

    def _get_result(self, md_file):
        """Format results as a dictionary."""
        log_file = open(md_file)
        lines = log_file.readlines()
        log_file.close()
        res = []
        append = False
        headline = None
        for l in lines:
            if append:
                res.append(l)
            if 'LJ-SR' in l:
                headline = l
                append = True
            if 'other-other' in l:
                break
        res_index = headline.split().index('LJ-SR') - 1
        res_dict = {}
        for l in res:
            res_dict[l.split()[0]] = l.split()[res_index]
        return res_dict


    def run(self, name, potential_type, verbose=False):
        """Compute RNAkb potential for a single file
        
        Args:

           path (str): name of a PDB file
           potential_type (str): '5pt' for 5 point or 'aa' for all atom aa is off

        Returns:

           a list of energies (strings) ['2.57116e+05', '1.62131e+03', '7.82459e+02', '3.00789e-01', '0.00000e+00', '0.00000e+00', '-2.51238e-03', '0.00000e+00', '2.59520e+05', '2.54174e-09', '2.59520e+05']

        .. warning:: 'aa' does not work because of "a line longer than 4095 characters"

        The function parses::

            Statistics over 1 steps using 1 frames
            Energies (kJ/mol)
            Bond          Angle    Proper Dih.  Improper Dih.          LJ-14
            2.44111e+05    1.76069e+03    8.12947e+02    1.82656e-01    0.00000e+00
            Coulomb-14        LJ (SR)   Coulomb (SR)      Potential    Kinetic En.
            0.00000e+00   -1.94624e-03    0.00000e+00    2.46685e+05    2.43227e-09
            Total Energy    Temperature Pressure (bar)
            2.46685e+05    6.67884e-10   -5.94232e+04
            Total Virial

            # ['Statistics', 'over', '1', 'steps', 'using', '1', 'frames', 'Energies', '(kJ/mol)', 
            'Bond', 'Angle', 'Proper', 'Dih.', 'Improper', 'Dih.', 'LJ-14', '2.44111e+05', '1.76069e+03', 
            '8.12947e+02', '1.82656e-01', '0.00000e+00', 'Coulomb-14', 'LJ', '(SR)', 'Coulomb', '(SR)',
            'Potential', 'Kinetic', 'En.', '0.00000e+00', '-1.94624e-03', '0.00000e+00', '2.46685e+05', 
            '2.43227e-09', 'Total', 'Energy', 'Temperature', 'Pressure', '(bar)', '2.46685e+05', '6.67884e-10', 
            '-5.94232e+04', 'Total', 'Virial']

            result[6] is really the RNAkb

        """
        prep = False

        self.sandbox_dir = self.sandbox_dir + os.sep

        self.log('Run %s' % name)
        if verbose: print('input: ', name, '\n', '-' * 80)
        copyfile(name, self.sandbox_dir + 'query.pdb')

        f = open(name)
        pdb_string = f.read()
        f.close()
        old_pdb = pdb_string
        pdb_string = make_rna_rnakb_ready(pdb_string)
        if old_pdb != pdb_string:
            self.pdb_fixes.append('make_rna_rnakb_ready')

        fn = self.sandbox_dir + 'query.pdb'
        if verbose: print('Gromacs ready file created: %s' % fn)
        f = open(fn, 'w')
        f.write(pdb_string)
        f.close()

        try:
            old_path = os.getcwd()
        except OSError:
            pass

        os.chdir(self.sandbox_dir)
        my_env = os.environ.copy()

        if potential_type == 'aa':
            my_env['GRMLIB'] = GRMLIB_PATH + ':' + DB_AA_PATH
            my_env['GMXLIB'] = GRMLIB_PATH + ':' + DB_AA_PATH
        elif potential_type == '5pt':
            my_env['GRMLIB'] = GRMLIB_PATH + ':' + DB_5PT_PATH
            my_env['GMXLIB'] = GRMLIB_PATH + ':' + DB_5PT_PATH
        else:
            raise Exception('RNAkb -- select aa or 5pt')

        if verbose:
            print('ff_aa:', GRMLIB_PATH + ':' + DB_AA_PATH)
            print('ff_5pt:', GRMLIB_PATH + ':' + DB_5PT_PATH)
            
        my_env['GMX_MAXBACKUP'] = '-1'
        my_env['LD_LIBRARY_PATH'] = GROMACS_LD_PATH

        # create topology files
        self.log('create topology files')
        # pdb2gmx
        out = run_command(self.executable[0], # '-v',
                       ['-quiet', '-ff', 'amber03', '-water', 'none' , '-f', self.sandbox_dir + 'query.pdb', '-o', self.sandbox_dir + 'query.gro',
                        '-p', self.sandbox_dir + 'query.top',
                        ],
                       env=my_env,
                       cmd_input='1\n6\n',
                       stdout_file=self.sandbox_dir + os.sep + 'out.txt',
                       stderr_file=self.sandbox_dir + os.sep + 'err.txt',
                       verbose=verbose,
                       )
        if out:
            print(out)
            self.log('Run failed')
            self.log_stdout_stderr()
            os.chdir(old_path)
            return -1

        # make index files
        fix_gromacs_gro(self.sandbox_dir + os.sep + 'query.gro')
        self.log('make index files')
        command = 'GRMLIB=' + my_env['GRMLIB'] + ' '
        command += 'GMXLIB=' + my_env['GMXLIB'] + ' '
        command += 'LD_LIBRARY_PATH=' + my_env['LD_LIBRARY_PATH'] + ' '


        # generate sandbox/groups.txt @groups
        gr = self.sandbox_dir + 'groups.txt'
        gtxt, energygrps, seq_uniq = prepare_groups(self.sandbox_dir + 'query.pdb', gr, potential=potential_type)
        if verbose: print('Groups file created: %s' % gr)

        # generate sandbox/score.mdp @score
        fout = self.sandbox_dir + 'score.mdp'

        if potential_type == 'aa': # ugly
            #format_score_mdp_aa(fout, energygrps, seq_uniq)
            format_score_mdp(fout, energygrps, seq_uniq)
        else:
            format_score_mdp(fout, energygrps, seq_uniq)

        command += self.executable[1] + ' -f query.gro > /dev/null 2> /dev/null < ' + gr
        #command += self.executable[1] + ' -f query.gro <' + gr
        bash_script = "#!/bin/bash\n" + command + "\n"
        with open(self.sandbox_dir + os.sep + 'make_index.sh', 'w') as f:            
            f.write(bash_script)
        cmd = 'bash ' + self.sandbox_dir + 'make_index.sh'# >/dev/null'
        if verbose: print('index cmd:', cmd)
        os.system(cmd)

        ## energy computation
        self.log('prepare energy computation')
        run_error = run_command(self.executable[3], # grompp
                       ['-noh', '-c', self.sandbox_dir + 'query.gro', '-p', self.sandbox_dir + 'query.top',
                        '-time', '1',
                        '-f', self.sandbox_dir + 'score.mdp',
                        '-n', self.sandbox_dir + 'index.ndx',
                        '-o', self.sandbox_dir + 'run.tpr', '-maxwarn', '15'],
                       env=my_env,
                       stdout_file=self.sandbox_dir + 'out.txt',
                       stderr_file=self.sandbox_dir + 'err.txt',
                       verbose=verbose,
                       )
        if not run_error:
            if verbose: print('prepare energy computation: OK', run_error)
        else:
            if verbose: print('prepare compute energy: FAIL, try to increase the box')
            # something went wrong, probably following lines will fix it
            f = open(self.sandbox_dir + os.sep + 'query.gro')
            box_size = f.readlines()[-1].split()
            f.close()

            box_size = [float(i) for i in box_size]
            max_size = 10

            for mult in range(5, max_size):
                ## increase box size
                if verbose: print('> Increase box size %d times' % mult)
                self.log('Increase box size %d times' % mult)

                run_error = run_command(self.executable[2], # editconf
                            ['-f', 'query.gro', '-o', 'query_big.gro',
                            '-box', ' '.join([str(i*mult) for i in box_size])],
                            env=my_env,
                            stdout_file=self.sandbox_dir + os.sep + 'out.txt',
                            stderr_file=self.sandbox_dir + os.sep + 'err.txt'
                            )

                if run_error:
                    self.log('Run failed', 'error')
                    self.log_stdout_stderr()
                    os.chdir(old_path)
                else:
                    if verbose: print('Increase box size: OK')

                # prepare energy computation again
                run_error =  run_command(self.executable[3], # grompp
                           ['-c', 'query_big.gro', '-p', 'query.top',
                            '-f', self.sandbox_dir + 'score.mdp',
                            '-n', self.sandbox_dir + 'index.ndx',
                            '-o', 'run.tpr', '-maxwarn', '15'],
                            env=my_env,
                            stdout_file=self.sandbox_dir + os.sep + 'out.txt',
                            stderr_file=self.sandbox_dir + os.sep + 'err.txt',
                            verbose=verbose)

                if run_error:
                    if verbose: print('Grompp: FAIL')
                    with open(self.sandbox_dir + os.sep + 'err.txt') as f:
                        txt = f.read()
                        if txt.find('referenced in the') > -1:
                            if verbose: print(txt) #print 'Group X referenced -- error'
                            return -2 # leave the loop
                        elif txt.find('The cut-off length is longer') > -1:
                            if verbose: print('The cuf-off length -- error')
                        else:
                            if verbose: print('unknown -- error, see ' + self.sandbox_dir + 'err.txt')
                            #sys.exit(1) # don't kill the function
                            pass

                    self.log_stdout_stderr()
                else:
                    if verbose: print('Grompp: OK')
                    break # leave the loop!

                if not os.path.exists(self.sandbox_dir + os.sep + 'run.tpr'):
                    self.log('Run failed %s' % name, 'error')
                    self.log_stdout_stderr()
                    os.chdir(old_path)

        self.log('compute energy')
        run_error = run_command(self.executable[4],
                       ['-table', 'table.xvg', '-tablep', 'tablep.xvg',
                        '-nt', '1',
                        '-s', self.sandbox_dir + 'run.tpr', '-e', self.sandbox_dir + 'run.edr'],
                       env=my_env,
                       stdout_file=self.sandbox_dir + os.sep + 'out.txt',
                       stderr_file=self.sandbox_dir + os.sep + 'err.txt',
                       verbose=verbose,
                       )
        if not run_error:
            if verbose: print('compute energy: OK', run_error)
        else:
            if verbose: print('compute energy: FAIL', run_error)
            self.log('Run failed', 'error')
            self.log_stdout_stderr()
            os.chdir(old_path)
            #return 255
            return [str(x) for x in [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]]
        try:
            os.chdir(old_pwd)
        except (OSError, NameError):
            pass

        # extract result and change format
        result = self._get_result(self.sandbox_dir + os.sep + 'md.log')
        result_sum = sum([float(result[k]) for k in result])

        md_txt = open(self.sandbox_dir + os.sep + 'md.log').read()
        energies_txt = re.search('Statistics over.*Total Virial', md_txt, flags=re.MULTILINE|re.DOTALL).group(0)
        energies = energies_txt.split()[16:21]+energies_txt.split()[29:34]+[energies_txt.split()[39]]
        if isnan(result_sum):
            self.log('Result is NaN', 'error')
            os.chdir(old_path)
            return [str(x) for x in [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]]
        else:
            self.log('Run successfull %s' % name)
            os.chdir(old_path)
            return energies
            #return result_sum != 0 and result_sum or -1

if __name__ == '__main__':
    wrapper = RNAkb(sandbox=True)

    fn = '1dqf_noCA.pdb' # problem, C is shifted to the right
    fn = '1dqf_noCA+ResnShift.pdb' # OK
    #fn = '1msy_output4_01-000001_AA.pdb' # OK
    #fn = 'cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_noC.pdb' # OK
    #fn = 'cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_error_fx.pdb'
    #fn = 'cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clx.pdb'
    #fn = 'C.pdb' 

    # 5pt
    if False:
        result = wrapper.run('test_data'
                             + os.sep + fn , '5pt', verbose=False) 
        print(fn, result)
        #wrapper.cleanup()

    # all
    result = wrapper.run('test_data' + os.sep + fn , '5pt', verbose=False) # or aa
    print(fn, result[6])
    #wrapper.cleanup()
