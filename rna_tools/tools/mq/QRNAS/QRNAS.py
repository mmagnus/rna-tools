#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing QRNAS.
https://github.com/sunandanmukherjee/QRNAS

Works both of M1 and Ubuntu.

Output::

    (...)

    Missing atom added: ATOM     34  H3T RA3 A  77      39.444  67.315  58.412  1.00  0.00           H
    Missing atom added: ATOM     23  OP3   G A   1      46.987  84.037  61.665  1.00  0.00
    Number of atoms:    2470
    Number of residues: 77
    Building intraresidual bonds...
    Building intraresidual angles...
    Building intraresidual dihedrals...
    Building intraresidual impropers...
    Building interresidual bonds...
    Building interresidual angles...
    Building interresidual dihedrals...
    Number of bonds built:     2662
    Number of angles built:    4764
    Number of dihedrals built: 7406
    Number of impropers built: 524
    -----------------------------------------------------------------------------
    End of molecule (chain)
    -----------------------------------------------------------------------------
    -----------------------------------------------------------------------------
    Building interresidual nonbonded pairs & H-bonds (it may take a while)...
    Number of electrostatic pairs:   0
    Number of van der Waals pairs:   332101
    Number of H-bonds built:         64
    Number of spring restraints:     0
    Number of positional restraints: 0
    Number of base pairs built:      20
    Number of molecules (chains) read: 1
    -----------------------------------------------------------------------------
    Performing minimization step: 1. Total energy = 28783.5671 kcal/mol (28753.2335 without restraints)
    Writing PDB file: /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpic_m3bcr/query_out.pdb ...Done.

"""

import os
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import QRNAS_PATH, QRNAS_CONFIG_PATH
import subprocess
import re
from shutil import copyfile

# directory where this script is
# files for some commands are also here
DIRECTORY = os.path.dirname(__file__)


class QRNAS(ProgramWrapper):
    """
    Wrapper class for AMBER force field via QRNA of Julisz Stasiewicz
    """
    program_name = 'QRNA'
    executable = 'QRNA'

    def __init__(self, job_id=None):
        super(QRNAS, self).__init__('sequence', 'seq_name', job_id=job_id)

    def _prepare_files(self):
        os.symlink(QRNAS_PATH + os.sep + self.executable,
                self.sandbox_dir + os.sep + 'QRNA')
        os.symlink(QRNAS_PATH + os.sep + 'forcefield',
                self.sandbox_dir + os.sep + 'forcefield')

    def run(self, filename, numSteps, verbose=False):
        """run()
        
        Return:
            e.g., [30307.1088', '28783.5671']

        """
        return [self.run_one(filename, str(numSteps), electrostatics=True, verbose=verbose), self.run_one(filename, str(numSteps), electrostatics=False, verbose=verbose)]

    def run_one(self, filename, numSteps, electrostatics=True, verbose=False):
        """
        Get AMBER energy via QRNA

        Args:
        path (str): The path of the file to wrap
        field_storage (FileStorage): The :class:Y instance to wrap
        temporary (bool): Whether or not to delete the file when the File
           instance is destructed

        Returns:
           BufferedFileStorage: A buffered writable file descriptor


        Arguments:
            * name = name (or path, if you're in a different directory)
              of a PDB file

        Output:
            * global energy as a floating point number
        """
        # copy file to sandbox and apply some fixes 
        copyfile(filename, self.sandbox_dir + os.sep + 'query.pdb')
        pdb_file = PDBFile(pdb_path=self.sandbox_dir + os.sep + 'query.pdb')
        #pdb_file.pedantic_pdb()
        #pdb_file.check_and_get_first_model()
        pdb_file.save(self.sandbox_dir + os.sep + 'query.pdb')
        self.pdb_fixes = pdb_file.fixes
        # write config
        f = open(self.sandbox_dir + os.sep + 'conf.cfg', 'w')
        f.write('NSTEPS ' + numSteps + '\nNUMTHREADS 1\n')
        if not electrostatics:
            f.write('ELECTR 0\n')
        f.close()

        #print self.sandbox_dir 
        cmd = self.sandbox_dir + os.sep + 'QRNA ' + \
                '-c  ' + self.sandbox_dir + os.sep + 'conf.cfg ' + \
                '-i ' + self.sandbox_dir + os.sep + 'query.pdb'
        if verbose:
            print(cmd)
        self.log(cmd, 'debug')
        self.log('Running program')
        out = subprocess.getoutput(cmd)
        self.log('Run finished')
        self.log(out, 'debug')

        if verbose: print(out)
        # get result
        #print out
        
        #x1 = re.compile('Number of electrostatic pairs:\s+(?P<energy>[0-9]+)').search(out).group('energy')
        #x2 = re.compile('Number of van der Waals pairs:\s+(?P<energy>[0-9]+)').search(out).group('energy')
        #x3 = re.compile('2Number of H-bonds built:\s+(?P<energy>[0-9]+)').search(out).group('energy')
        
        rx = re.compile('Performing minimization step: ' + numSteps + '. Total energy = (?P<energy>[01-9\.e\+]+) kcal').search(out)
        if rx:
            energy = rx.group('energy')
        else:
            energy = -1
            self.log('# problem --', 'error')
        return energy


def main():
    qrnas = QRNAS()
    try:
        # Performing minimization step: 1. Total energy = 28783.5671 kcal/mol (28753.2335 without restraints)
        filename = "test" + os.sep + "unmod_Val3_tRNA_model_si.pdb"
        energy = qrnas.run(filename, 1, False)
        print('energy: ', energy)
    except Exception as e:
        print(e)
    finally:
        #qrna.cleanup()
        pass

if __name__ == "__main__":
    main()
