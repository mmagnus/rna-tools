#!/usr/bin/python

import subprocess
import sys
import os


def test_cmd():
    curr = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cmd = "./rna_pdb_merge_into_one.py test_in/*.pdb > test_out/rp17.pdb"
    o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    assert out == ''
    assert err == ''
    os.chdir(curr)


test_cmd()
