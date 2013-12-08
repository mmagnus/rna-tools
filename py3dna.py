#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Parser to 3dna <http://x3dna.org/>

Installation::

  BINARY_PATH = '/usr/bin/x3dna-dssr-64bit'

Usage::

    $ python py3dna.py
    test_data/1xjr.pdb
    gGAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU
    ..(((((((...((((.((((.....))..))..))).).)))))))
    test_data/6TNA.pdb
    GCGGAUUUAgCUCAGuuGGGAGAGCgCCAGAcUgAAgAPcUGGAGgUCcUGUGuPCGaUCCACAGAAUUCGCACCA
    (((((((..((((.....[..)))).((((.........)))).....(((((..]....))))))))))))....

.. warning:: it get_seq() and get_secstruc work only for 1 chain PDB file!

"""

from subprocess import Popen, PIPE
from os import listdir, remove

import re

BINARY_PATH = '/usr/bin/x3dna-dssr-64bit'
BINARY_PATH_FP = '/home/magnus/opt/x3dna-dssr/2.1/x3dna-v2.1/bin/find_pair'


class Py3DNAMissingFile(Exception):

    pass


class Py3DNA(object):

    """
    Atributes:

     **curr_fn**
     **report**

    """

    def __init__(self, pdbfn):
        """Set self.curr_fn based on pdbfn"""

        self.curr_fn = pdbfn
        self.run_3dna()
        self.clean_up()

    def get_report(self):
        """Run find_pair

        ..warning:: To get modification use get_modifications()
        """
        cmd = BINARY_PATH_FP + ' ' + self.curr_fn + ' stdout | analyze stdin'  # hack
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = out.stdout.read()
        outerr = out.stderr.read()

        text = ''
        for l in outerr.split('\n'):
            if l.startswith('Time used') or l.startswith('..') or l.startswith('handling') or l.startswith('uncommon residue'):
                continue
            if l.strip():
                text += l + '\n'
        return text.strip()

    def get_modifications(self):
        """Run find_pair to find modifications.
        """
        cmd = BINARY_PATH_FP + ' -p ' + self.curr_fn + ' /tmp/fpout'  # hack
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        outerr = out.stderr.read()

        text = ''
        for l in outerr.split('\n'):
            if l.startswith('uncommon residue'):
                text += l.replace('uncommon residue ', '') + '\n'
        return text.strip()

    def run_3dna(self):
        """
        """

        cmd = BINARY_PATH + ' -i=' + self.curr_fn
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = out.stdout.read()
        outerr = out.stderr.read()

        f = open('py3dna.log', 'w')
        f.write(cmd + '\n' + stdout)
        f.close()

        if outerr.find('does not exist!') > -1:  # not very pretty
            raise Py3DNAMissingFile
        if outerr.find('not found') > -1:  # not very pretty
            raise Exception('x3dna not found!')

        rx = re.compile('no. of DNA/RNA chains:\s+(?P<no_DNARNAchains>\d+)\s').search(stdout)
        if rx:
            no_of_DNARNA_chains = int(rx.group('no_DNARNAchains'))
            msg = 'py3dna::no_of_DNARNA_chains'
            self.report = msg + '\n' + msg + '\n' # hack!
        else:
            raise Exception('no_of_DNARNA_chains not found')

        if no_of_DNARNA_chains:
            self.report = stdout.strip()

    def get_ion_water_report(self):
        """@todo
File name: /tmp/tmp0pdNHS
    no. of DNA/RNA chains: 0 []
    no. of nucleotides:    174
    no. of waters:         793
    no. of metals:         33 [Na=29, Mg=1, K=3]
        """
        pass

    def clean_up(self, verbose=False):
        files_to_remove = [
            'dssr-helices.pdb',
            'dssr-pairs.pdb',
            'dssr-torsions.dat',
            'dssr-hairpins.pdb',
            'dssr-multiplets.pdb',
            'dssr-stems.pdb',
            'dssr-Aminors.pdb',
            ]

        for f in files_to_remove:
            try:
                remove(f)
            except OSError:
                if verbose: print 'can not remove %s' % f

    def get_seq(self):
        """Get sequence.

        Somehow 1bzt_1 x3dna	UCAGACUUUUAAPCUGA, what is P?
        P -> u
        """

        return self.report.split('\n')[-2].replace('P','u').replace('I','a')

    def get_secstruc(self):
        """Get secondary structure.
        """

        return self.report.split('\n')[-1]


if __name__ == '__main__':
    directory = 'test_data'

    pdbs = listdir(directory)
    #pdbs = ['1xjr.pdb']
    for f in pdbs:
        fn = directory + '/' + f
        print fn

        p = Py3DNA(fn)

        s = p.get_seq()
        print s

        s = p.get_secstruc()
        print s
