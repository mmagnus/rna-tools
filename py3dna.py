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

BINARY_PATH = '/usr/bin/x3dna-dssr-64bit'


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

    def run_3dna(self):
        """
        """

        cmd = BINARY_PATH + ' -i=' + self.curr_fn
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = out.stdout.read()
        outerr = out.stderr.read()

        if outerr.find('does not exist!') > -1:  # not very pretty
            raise Py3DNAMissingFile
        if outerr.find('not found') > -1:  # not very pretty
            raise Exception('x3dna not found!')

        self.report = stdout.strip()

    def clean_up(self):
        files_to_remove = [
            'dssr-helices.pdb',
            'dssr-pairs.pdb',
            'dssr-torsions.dat',
            'dssr-hairpins.pdb',
            'dssr-multiplets.pdb',
            'dssr-stems.pdb',
            ]

        for f in files_to_remove:
            remove(f)

    def get_seq(self):
        """Get sequence.
        """

        return self.report.split('\n')[-2]

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
