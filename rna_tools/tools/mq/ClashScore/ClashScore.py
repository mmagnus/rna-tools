#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
ic.configureOutput(prefix='> ')

import argparse
import os

from rna_tools.rna_tools_config import PHENIX_BIN_PATH
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from subprocess import Popen, PIPE

# directory where this script is
# files for some commands are also here
DIRECTORY = os.path.dirname(__file__)


class ClashScore(object):
    """ClashScore: Wrapper class for running phenix.clashscore
    """
    def __init__(self):
        pass

    def run(self, fn, verbose=False):
        """
        Args:

            fn (string): path to a file
            verbose (bool): be verbose

        Returns:
             
            clashscore (float)

        Example::

            /Applications/phenix-1.18.2-3874/build/bin/phenix.clashscore test/1xjrA.pdb
            Using electron cloud x-H distances and vdW radii

            Adding H/D atoms with reduce...

            Bad Clashes >= 0.4 Angstrom:
             A  17    A  N1   A  34    G  N1  :0.430
            clashscore = 0.66
            test/1xjrA.pdb
            0.66
            /Applications/phenix-1.18.2-3874/build/bin/phenix.clashscore test/1xjrA_M1.pdb
            Using electron cloud x-H distances and vdW radii

            Adding H/D atoms with reduce...

            Bad Clashes >= 0.4 Angstrom:
             A  41    G  H2'  A  42    U  OP1 :0.738
             A  25    U  H6   A  25    U HO2' :0.659
             A  46    U  H2'  A  47    U  OP2 :0.656
             A  29    A H5''  A  30    U  O2' :0.623
             (...)

              A  24    G  O2'  A  26    A  C8  :0.410
             A  42    U  O2'  A  43    G  H8  :0.409
             A  43    G  H2'  A  44    A  O4' :0.408
             A  28    G  C8   A  28    G  O2' :0.403
            clashscore = 15.45
            test/1xjrA_M1.pdb
            15.45
            
            # so the output is
            clashscore = 15.45

        """
        keep_hydrogens = False
        if verbose: ic(keep_hydrogens)
        if keep_hydrogens:
            cmd = PHENIX_BIN_PATH + os.sep + 'phenix.clashscore %s keep_hydrogens=True ' % fn
        else:
            cmd = PHENIX_BIN_PATH + os.sep + 'phenix.clashscore %s ' % fn
        if verbose: print(cmd)
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)
        stdout = out.stdout.read().decode().strip()
        if verbose: print(stdout)
        out = float(stdout.split('\n')[-1].replace('clashscore = ',''))

        # this is another 
        #cmd = 'phenix.clashscore_1.8.1-1168 keep_hydrogens=True %s' % name
        #out2 = float(getoutput(cmd).split('\n')[-1].replace('clashscore = ',''))
        
        #return ','.join([str(out), str(out2)]) 
        return out

    def cleanup(self):
        pass
    

def test():
    wrapper = ClashScore()
    fns = ['../test' + os.sep + '1xjrA.pdb', 
           '../test' + os.sep + '1xjrA_M1.pdb',  # native, and M1 from rasp decoys
           '../test' + os.sep + '6TNA.pdb'
    ]
    for f in fns:
        result = wrapper.run(f, True)
        print(f)
        print(result)
        wrapper.cleanup()

if __name__ == '__main__':
    test()
