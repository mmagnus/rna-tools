#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing Dfire potential

Installation:

     git clone https://github.com/tcgriffith/dfire_rna.git
     make
     # add DFIRE_RNA_HOME to .bashrc

"""
from __future__ import print_function
import argparse
import os
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import dfire_PATH

class Dfire(ProgramWrapper):
    """
    Wrapper class for Dfire.
    """
    def __init__(self):
        super(Dfire, self).__init__()

    def run(self, path_to_pdb, verbose=False):
        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        old_pwd = os.getcwd()
        os.chdir(self.sandbox_dir)
        self.log('dfire::start for %s' % self.sandbox_dir + '/query.pdb')
        cmd = "export DFIRE_RNA_HOME=" + dfire_PATH + "; "
        cmd += dfire_PATH + '/bin/DFIRE_RNA ' + self.sandbox_dir + '/query.pdb >'  + self.sandbox_dir + '/log.txt.dfire.rna'
        if verbose: print(cmd)
        os.system(cmd)

        self.log('dfire::Run finished')

        for line in open(self.sandbox_dir + '/log.txt.dfire.rna'):
            score = line.strip().split()[1] # /tmp/tmplDNztf/query.pdb -12480.188898
        if verbose: print(open(self.sandbox_dir + '/log.txt.dfire.rna').read())
        os.chdir(old_pwd)
        return float(score)


def get_parser():
        parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

        parser.add_argument("-v", "--verbose",
                            action="store_true", help="be verbose")
        parser.add_argument("--test", action="store_true")
        parser.add_argument("--file", help="", default="")  # nargs='+')
        return parser


def test(verbose):
    wrapper = Dfire()
    try:
        result = wrapper.run('../test' + os.sep + '1a9n.pdb', verbose)
        if result:
            print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

    print((wrapper.run('../../../input/mq/' + os.sep + '1a9nR.pdb', verbose)))
    print((wrapper.run('../../../input/mq/' + os.sep + '3b58ABC.pdb', verbose)))
    print((wrapper.run('../../../input/mq/' + os.sep + '5e3hBC.pdb', verbose)))
    print((wrapper.run('../../../input/mq/' + os.sep + 'S_000001_000.pdb', verbose)))
    

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
        
    if args.test:
        test(args.verbose)
    else:
        wrapper = Dfire()
        result = wrapper.run(args.file, args.verbose)
        print(result)
