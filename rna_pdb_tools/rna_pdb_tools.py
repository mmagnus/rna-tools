#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""rna_pdb_tools main parser code."""

import argparse
import os
import time

from pdb_parser_lib import *

def get_parser():
    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()
    parser = argparse.ArgumentParser('rna-pdb-tools.py ver: %s' % version)

    parser.add_argument('-r', '--report', help='get report',
                        action='store_true')

    parser.add_argument('-c', '--clean', help='get clean structure',
                        action='store_true')

    parser.add_argument('--get_chain', help='get chain, .e.g A')

    parser.add_argument('--fetch', action='store_true', help='fetch file from the PDB db')

    parser.add_argument('--fetch_ba', action='store_true', help='fetch biological assembly from the PDB db')

    parser.add_argument('--get_seq', help='get seq', action='store_true')

    parser.add_argument('--rosetta2generic', help='convert ROSETTA-like format to a generic pdb',
                        action='store_true')

    parser.add_argument('--get_rnapuzzle_ready', help='get RNApuzzle ready (keep only standard atoms, renumber residues)',
                        action='store_true')

    parser.add_argument('--no_hr', help='do not insert the header into files',
                        action='store_true')

    parser.add_argument('--renumber_residues', help='',
                        action='store_true')

    parser.add_argument('--collapsed_view', help='',
                        action='store_true')

    parser.add_argument('--cv', help='alias to collapsed_view',
                        action='store_true')

    parser.add_argument('--edit',
			dest="edit",
                        default='',
                        help="edit 'A:6>B:200', 'A:2-7>B:2-7'")

    parser.add_argument('--delete',# type="string",
			dest="delete",
			default='',
			help="delete the selected fragment, e.g. A:10-16")
    
    parser.add_argument('file', help='file')
    #parser.add_argument('outfile', help='outfile')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if args.report:
        s = StrucFile(args.file)
        print(s.get_report())
        print(s.get_preview())
        print(s.get_info_chains())

    if args.clean:
        s = StrucFile(args.file)
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        #print s.get_preview()
        #s.write(args.outfile)
        if not args.no_hr:
            add_header()
        print(s.get_text())

    if args.get_seq:
        s = StrucFile(args.file)
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        #print s.get_preview()
        print(s.get_seq())

    if args.get_chain:
        s = StrucFile(args.file)
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        #print s.get_preview()
        print(s.get_chain(args.get_chain))

    if args.rosetta2generic:
        s = StrucFile(args.file)
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.fix_op_atoms()
        s.renum_atoms()
        #print s.get_preview()
        #s.write(args.outfile)
        if not args.no_hr:
            add_header()
        print(s.get_text())

    if args.get_rnapuzzle_ready:
        s = StrucFile(args.file)
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.fix_op_atoms()
        s.renum_atoms()
        #print s.get_preview()
        #s.write(args.outfile)
        if not args.no_hr:
            add_header()
        s.get_rnapuzzle_ready(args.renumber_residues)
        print(s.get_text())

    if args.renumber_residues:
        s = StrucFile(args.file)
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.get_rnapuzzle_ready(args.renumber_residues)
        s.renum_atoms()
        if not args.no_hr:
            add_header()
        print(s.get_text())

    if args.delete:
        selection = select_pdb_fragment(args.delete)
        s = StrucFile(args.file)
        if not args.no_hr:
            add_header()
            print('HEADER --delete ' + args.delete) #' '.join(str(selection))
        for l in s.lines:
            if l.startswith('ATOM'):
                chain = l[21]
                resi = int(l[23:26].strip())
                if chain in selection:
                    if resi in selection[chain]:
                        continue  # print chain, resi
                print(l)

    if args.edit:
        edit_pdb(args)
        
    if args.fetch:
        fetch(args.file)

    if args.fetch_ba:
        fetch_ba(args.file)

    if args.collapsed_view or args.cv:
        collapsed_view(args)

if __name__ == '__main__':
    pass
