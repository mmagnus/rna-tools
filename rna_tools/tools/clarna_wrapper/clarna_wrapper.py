#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Wrapper for ClaRNA <http://genesilico.pl/clarna/>
"""

import tempfile
import os
import sys

from .clarna_wrapper_config import PATH
from subprocess import Popen, PIPE

VERBOSE = False
DEV = False

class ClaRNAWrapper(object):
    def __init__(self, pdbfn):
        """Set self.curr_fn based on pdbfn
        """
        self.name = os.path.basename(pdbfn).replace('.pdb', '')
        self.curr_fn = pdbfn
        self.run_clarna()
        self.run_gts()
        if not DEV: self.clean_up()

    def run_clarna(self):
        """
        """
        self.tmp = tempfile.NamedTemporaryFile().name + '.clarna'

        if not os.path.isfile(os.path.abspath(self.curr_fn)):
            raise Exception('Error: File does not exist!')

        cmd = PATH + '/clarna.py --save-graph="' + self.tmp + '" ' + self.curr_fn
        if VERBOSE: print('exe: %s' % cmd)
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = out.stdout.read()
        outerr = out.stderr.read()

        self.log = self.tmp + '.clarna.log'
        if VERBOSE: print('log:', self.log)
        f = open(self.log, 'w')
        f.write(cmd + '\n' + stdout)
        f.close()

        self.report = stdout.strip()

    def run_gts(self):
        # run  graph-to-secondary-structure.py
        self.tmp_ss = self.tmp + '.ss'

        cmd = PATH + '/graph-to-secondary-structure.py -i "' + self.tmp + '" -o "' + self.tmp_ss + '"'
        if VERBOSE: print('exe: %s' % cmd)
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = out.stdout.read()
        outerr = out.stderr.read()
        if outerr:
            raise Exception(outerr)

        self.log_ss = self.tmp + '.clarna.ss.log'
        if VERBOSE: print('log:', self.log_ss)
        f = open(self.log_ss, 'w')
        f.write(cmd + '\n' + stdout)
        f.close()

        self.gts = open(self.tmp_ss).read().strip()
        if VERBOSE: print(self.gts)

    def clean_up(self, verbose=False):
        files_to_remove = [
            self.tmp,
            self.log,
            self.tmp_ss,
            self.log_ss,
            ]

        for f in files_to_remove:
            try:
                os.remove(f)
            except OSError:
                raise Exception('File can not be remove %s' % f)

    def get_seq(self):
        return self.gts.split()[-2]

    def get_ss(self):
        return self.gts.split()[-1]

if __name__ == '__main__':
    sys.argv.append('')
    sys.argv[1] = 'test_data/receptor.pdb'
    try:  ## quick and dirty
        f = sys.argv[1]
    except IndexError:
        print('clarna_wrapper.py <fn>')
        sys.exit(1)

    c = ClaRNAWrapper(f)
    print('>', c.name)
    print(c.get_seq())
    print(c.get_ss())
