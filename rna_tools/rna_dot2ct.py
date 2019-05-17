#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The output file is <input-file>.ct

Wrapper to

RNAstructure: software for RNA secondary structure prediction and analysis. (2010). RNAstructure: software for RNA secondary structure prediction and analysis., 11, 129. http://doi.org/10.1186/1471-2105-11-129
"""
from __future__ import print_function
import textwrap
import argparse
import subprocess

def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err

def dot2ct():
    pass


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help=textwrap.dedent("""Input is:
>seq
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
((...((((((((((.......)))))))))).))
"""))
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    cmd='dot2ct ' + args.file + ' ' + args.file + '.ct'
    out, err= exe(cmd)
    print(out, end='')
    if not err:
        print(' Created: %s' % args.file + '.ct')
