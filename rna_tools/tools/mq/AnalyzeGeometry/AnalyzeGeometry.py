#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing AnalyzeGeomotery.
"""
import os
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile  #resname_check_and_3to1, set_residues_bfactor
from rna_tools.rna_tools_config import PHENIX_BIN_PATH
from subprocess import Popen, PIPE

# directory where this script is
# files for some commands are also here
DIRECTORY = os.path.dirname(__file__)


class AnalyzeGeometry(object):
    """
    Wrapper class for running clashscore

    .. note::  Sequence is required to get the length to calculate % of corrent residues.
    """
    def __init__(self, verbose=False):
        pass

    def run(self, name, verbose=False):
        """
        Args:

           name (str): the path of the file to wrap
           verbose (boolen): be verbose

        Returns:

            float: score; % of Backbone torsion suites (# of them per seq)

        Output::

            ----------Backbone torsion suites----------

              Suite ID                 suite  suiteness  triaged angle
                 A A   3                  !!      0.000     delta
                 G A   4                  !!      0.000  epsilon-1
                 G A  11                  !!      0.000      None
                 C A  20                  !!      0.000     gamma
                 U A  25                  !!      0.000     delta
                 A A  26                  !!      0.000   delta-1
                 A A  29                  !!      0.000      None
                 U A  30                  !!      0.000  epsilon-1
                 G A  41                  !!      0.000     delta
                 U A  42                  !!      0.000   delta-1
                 A A  45                  !!      0.000     delta
                 U A  46                  !!      0.000  epsilon-1
                 U A  47                  !!      0.000   delta-1
              11 suites triaged and 0 incomplete leaving 35 suites
              13/46 suite outliers present
              Average suiteness: 0.490
            13 # count lines after 'traged angle' and minus 3 (# of last lines)
            test/1xjrA_M1.pdb
            34.7826

        Output for a perfect structure::

            /Applications/phenix-1.18.2-3874/build/bin/phenix.rna_validate /Users/magnus/work/src/rna-tools/rna_tools/input/mq/5e3hBC.pdb
            CGACGCUAGCGUACGCUAGCGUCG
            AnalyzeGeomtery:: ----------Backbone bond lenths----------                   

              All bonds within 4.0 sigma of ideal values.

                                ----------Backbone bond angles----------                   

              All angles within 4.0 sigma of ideal values.

                                    ----------Sugar pucker----------                       

              All puckers have reasonable geometry.

                              ----------Backbone torsion suites----------                  

              0 suites triaged and 0 incomplete leaving 24 suites
              All RNA torsion suites are reasonable.
              Average suiteness: 0.766

        """ 
        #cmd = 'phenix.rna_validate_1.8.1-1168 %s' % name
        cmd = PHENIX_BIN_PATH + os.sep + 'phenix.rna_validate %s' % name

        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)
        out = out.stdout.read().decode().strip()

        c = 0
        now_count = False

        pdb_file = PDBFile(pdb_path=name)
        seq = pdb_file.seq_from_pdb()
        if verbose:
            print(cmd)
            print(seq)
            print('AnalyzeGeomtery::', out)

        for l in out.split('\n'):
            if now_count:
                c += 1
            if 'triaged angle' in l: # #suiteID:suite:suiteness:triaged_angle'):
                now_count = True
        # if c = 0 it means that triagged angle was not even detected
        if c:
            c = c - 3
        return round(c/float(len(seq))*100,4)
        #return ','.join([str(out), str(out2)]) 

    def cleanup(self):
        pass


def main():
    wrapper = AnalyzeGeometry()
    fns = [#'test' + os.sep + '1xjrA.pdb', # <- gtp causes error
           '../test' + os.sep + '1xjrA_M1.pdb',
           #'test' + os.sep + '1xjrA_M500.pdb'
    ] # native, and M1 from rasp decoys
    for f in fns:
        result = wrapper.run(f, True)
        print(f)
        print(result)
        wrapper.cleanup()

if __name__ == '__main__':
    main()
