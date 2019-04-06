#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for rebuilding all atom representation from reduced,
3 atom representation.
"""
import os
DIRNAME = os.path.dirname(__file__)

from subprocess import Popen
from shutil import copyfile, rmtree
from tempfile import mkdtemp
from .SingleLineUtils import get_res_num
from .PDBFile import PDBFile

REBUILDRNA_PATH=os.path.join(DIRNAME, 'RebuildRNA', 'trunk')

def check_3_atom(pdb_string):
    '''Check what number of atoms per residue in PDB

    Arguments:
      * pdb_string = string with PDB structure to check

    Output:
      * True if there are 3 atoms per residue, False otherwise
    '''
    pdb_lines = pdb_string.split('\n')
    atom_lines = [l for l in pdb_lines if l.startswith('ATOM')]
    res_nums = [get_res_num(l) for l in atom_lines]
    uniq_nums = list(set(res_nums))
    counted = [res_nums.count(i) for i in uniq_nums]
    return not [i for i in counted if i != 3]


def rebuild_full_atom(input_pdb, output_pdb):
    """Make all atom PDB file from 3 atom PDB file

    Arguments:
      * input_pdb = path to input PDB
      * output_pdb = path to output PDB
    """
    try:
        old_pwd = os.getcwd()
        print(old_pwd)
    except:
        old_pwd = None
    tempdir = mkdtemp()
    os.chdir(tempdir)
    try:
        os.symlink(REBUILDRNA_PATH + os.sep + 'histograms', 'histograms')
        os.symlink(REBUILDRNA_PATH + os.sep + 'pdb', 'pdb')
    except OSError:
        pass
    # remove strange residue names
    pdb_file = PDBFile(pdb_path=os.path.join(old_pwd, input_pdb))
    pdb_file._resname_3to1()
    pdb_file.save(os.path.join(tempdir, input_pdb))
    cmd = Popen([REBUILDRNA_PATH + os.sep + 'RebuildRNA',
        '-t', os.path.join(tempdir, input_pdb), '-o', tempdir])
    returncode = cmd.wait()
    out_name = tempdir + os.sep +\
        '.'.join(os.path.basename(input_pdb).split('.')[:-1]) + '_beforeCorrecting.pdb'
    os.rename(out_name, output_pdb)
    try:
        if old_pwd:
            os.chdir(old_pwd)
        rmtree(tempdir)
    except OSError:
        pass
    return returncode


def check_and_rebuild(input_pdb, output_pdb):
    """Check PDB representation and rebuild if necessary

    Arguments:
      * input_pdb = path to input PDB
      * output_pdb = path to output PDB

    Output:
      * True if rebuilding was necessary, False if not
    """
    f = open(input_pdb)
    input_pdb_string = f.read()
    f.close()
    if check_3_atom(input_pdb_string):
        rebuild_full_atom(input_pdb, output_pdb)
        return True
    else:
        copyfile(input_pdb, output_pdb)
        return False
