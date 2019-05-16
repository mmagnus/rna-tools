#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
output is <file>.ct
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
