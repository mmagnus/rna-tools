#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
ATOM    561  O3'   C    27
"""
import sys
import os
import argparse

PATH = os.path.abspath(__file__)
if os.path.islink(PATH):
    PATH = os.path.dirname(os.readlink(PATH)) + '/' # / don't forget
else:
    PATH = os.path.dirname(os.path.abspath(__file__)) + '/'

try:
    from diffpdb_conf import DIFF_TOOL
except ImportError:
    DIFF_TOOL = 'diff' ## kompare

def do_file_atom_res(fn):
    """Take a file path, open the file, for each line take 
    first 31 characters, saves these lines to a new file."""

    text_new = ''
    for l in open(fn):
        if l.startswith('ATOM') or l.startswith('HETATM') or l.startswith('END') or l.startswith('TER') or l.startswith('MODEL') or l.startswith('ENDMDL'):
            text_new +=  l[12:20].strip() + '\n'
    open(fn + '.out', 'w').write(text_new)

def do_file_atom_res_and_resi(fn):
    """Take a file path, open the file, for each line take 
    first 31 characters, saves these lines to a new file."""

    text_new = ''
    for l in open(fn):
        if l.startswith('ATOM') or l.startswith('HETATM') or l.startswith('END') or l.startswith('TER') or l.startswith('MODEL') or l.startswith('ENDMDL'):
            text_new +=  l[12:26].strip() + '\n'
    open(fn + '.out', 'w').write(text_new)
#ATOM   1002  P     A A  48

def do_file(fn):
    """Take a file path, open the file, for each line take 
    first 31 characters, saves these lines to a new file."""

    text_new = ''
    for l in open(fn):
        if l.startswith('ATOM') or l.startswith('HETATM') or l.startswith('END') or l.startswith('TER') or l.startswith('MODEL') or l.startswith('ENDMDL'):
            l = l[:31].strip()
            text_new += l + '\n'
    open(fn + '.out', 'w').write(text_new)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('diffpdb.py')
    parser.add_argument('--names', help='take only atom residues names', action='store_true')
    parser.add_argument('--names_and_resi', help='take only atom residues names', action='store_true')
    parser.add_argument('--htmlout', help='take only atom residues names', action='store_true')
    parser.add_argument('--method', help='method e.g. `diff`')#, action='store_true')
    parser.add_argument('f1', help='file')     
    parser.add_argument('f2', help='file')     

    args = parser.parse_args()
    
    str1_fn = args.f1
    str2_fn = args.f2

    if args.method:
        DIFF_TOOL = args.method

    if args.names:
        do_file_atom_res(str1_fn)
        do_file_atom_res(str2_fn)
    elif args.names_and_resi:
        do_file_atom_res_and_resi(str1_fn)
        do_file_atom_res_and_resi(str2_fn)
    else:
        do_file(str1_fn)
        do_file(str2_fn)

    if args.htmlout:
        cmd = 'diff -u %s %s | %s/lib/diff2html.py' % (str1_fn + '.out', str2_fn + '.out', PATH)
        os.system(cmd)
    else:
        os.system(' '.join([DIFF_TOOL, str1_fn + '.out', str2_fn + '.out']))
