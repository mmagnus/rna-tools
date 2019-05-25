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
import tempfile
import six


def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err

def dot2ct_file(path, path_output='', verbose=False):
    """
    Args:
        path (str): a path to the input file
        path_output (str): a path to the output file, by default it's '<path>' + '.ct'

    Returns:
        path to the output file
    """
    if not path_output:
        path_output = path + '.ct'
    cmd='dot2ct ' + path + ' ' + path_output + '.ct'
    out, err= exe(cmd)
    if verbose: print(out, end='')
    if not err:
        if verbose: print(' Created: %s' % args.file + '.ct')
    return path_output


def dot2ct(seq, ss, verbose=True):
    """
    Args:
        seq (str|RNAstructure object): sequence
        ss  (str): secondary structure

    Returns:
        ct (str): content of the output ct file
    """
    t = tempfile.NamedTemporaryFile(delete=False)
    with open(t.name, 'w') as f:
        f.write('> ss\n')
        if not isinstance(seq, six.string_types):
            seq = seq.seq
        f.write(seq.strip() + '\n')
        f.write(ss.strip() + '\n')

    cmd='dot2ct ' + t.name + ' ' + t.name + '.ct'
    if verbose: print(cmd)
    out, err = exe(cmd)
    print(out, end='')
    if not err:
        print(' Created: %s' % f + '.ct')
    else:
        print(err)
    with open(t.name + '.ct') as f:
        return f.read().strip()


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

    dot2ct_file(args.file)
