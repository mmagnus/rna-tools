#!/usr/bin/python
# -*- coding: utf-8 -*-

"""diffpdb - a simple tool to compare text-content of PDB files

The method is quick-and-dirty, but works!

The script takes first 31 characters of lines (or only atom names and residue names)
starting with ``HETATM`` or ``ATOM`` and save these lines to a <filename>.out file.

One file is created per pdb. In the final step DIFF_TOOL is executed
on these two output files. You get a diff output. That's it! Enjoy!

Configuration:

 * ``DIFF_TOOL="open -a diffmerge"`` or ``DIFF_TOOL="kompare"`` to set up what tool would you like to use to diff files in the file ``rna-pdb-tools/tools/diffpdb/diffpdb_conf.py`` (create it if needed)

.. image:: ../../rna_tools/tools/diffpdb/doc/screenshot.png

``./diffpdb.py --names test_data/4/1duq.pdb test_data/4/1duq_decoy0171_amb_clx.pdb``

.. image:: ../../rna_tools/tools/diffpdb/doc/screenshot2.png

and on the Mac (using ``diffmerge``):

.. image:: ../../rna_tools/tools/diffpdb/doc/diffpdb_osx_diffmerge.png

One of the difference that can be detected with the script is variants of atoms.

.. image:: ../../rna_tools/tools/diffpdb/doc/atom-variants.png

or a detection of missing atom.

.. image:: ../../rna_tools/tools/diffpdb/doc/missing-atoms.png

or a detection of missing OP3 at the beginning.

.. image:: ../../rna_tools/tools/diffpdb/doc/missing-op3.png

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
    from rna_tools.rna_tools_config import DIFF_TOOL
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

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__ + '\n' + '  method: %s ' % DIFF_TOOL, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--names', help='take only atom residues names', action='store_true')
    parser.add_argument('--names-and-resi', help='take only atom residues names', action='store_true')
    parser.add_argument('--htmlout', help='take only atom residues names', action='store_true')
    parser.add_argument('--method', help='method e.g. `diff`')#, action='store_true')
    parser.add_argument('--verbose', help='be verbose', action='store_true')
    parser.add_argument('f1', help='file')
    parser.add_argument('f2', help='file')
    return parser

if __name__ == '__main__':
    parser = get_parser()
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
        # -s
        #cmd = 'diff -u %s %s | %s/lib/diff2html.py' % (str1_fn + '.out', str2_fn + '.out', PATH)
        cmd = 'diff %s %s ' % (str1_fn + '.out', str2_fn + '.out', PATH)
        os.system(cmd)
    else:
        # DIFF_TOOL
        cmd = ' ' + ' '.join(['diff', ' -s -y ', str1_fn + '.out', str2_fn + '.out'])
        if args.verbose:
            print(cmd)
        os.system(cmd)
