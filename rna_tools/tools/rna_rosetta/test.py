#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import sys


def test():
    cmd = "./rna_rosetta_run.py -r -g -n 10 -c 1 --sandbox test_data test_data/tetraloop.fa"
    # - i - e
    o = subprocess.Popen(cmd, shell=True, stdout=sys.stdout.flush(), stderr=subprocess.PIPE)
    stderr = o.stderr.read().strip()
    assert stderr == ''

def test_cluster():
    cmd = "./rna_rosetta_cluster.py test_data/selected.out 10"
    o = subprocess.Popen(cmd, shell=True, stdout=sys.stdout.flush(), stderr=subprocess.PIPE)
    stderr = o.stderr.read().strip()
    assert stderr == ''



if __name__ == '__main__':
    test()
    test_cluster()
