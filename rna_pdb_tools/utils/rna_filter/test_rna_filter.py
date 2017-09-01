#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import subprocess
import git

def run_cmd(cmd, verbose=True):
    if verbose:
        print(cmd)
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err

def is_modifed(fn):
    out, err = run_cmd('git status %s' % fn)
    if 'modified' in out or 'untracked file' in out:
        return True
    else:
        return False

def test_rna_filter():
    fn = "test_data/rna_filter.txt"
    cmd = "python rna_filter.py -r test_data/restraints.txt -s test_data/CG.pdb | tee test_data/rna_filter.txt"
    out, err = run_cmd(cmd)
    assert is_modifed(fn) == False
