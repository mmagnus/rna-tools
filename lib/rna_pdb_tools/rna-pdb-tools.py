#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import time

from pdb_parser_lib import *
from rna_pdb_tools.utils.extra_functions.select_fragment import select_pdb_fragment_pymol_style, select_pdb_fragment

def add_header():
    now = time.strftime("%c")
    print 'HEADER Generated with rna-pdb-tools'
    print 'HEADER ver %s \nHEADER https://github.com/mmagnus/rna-pdb-tools \nHEADER %s' % (version, now)

if __name__ == '__main__':
    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()
    parser = argparse.ArgumentParser('rna-pdb-tools.py ver: %s' % version)

    parser.add_argument('-r', '--report', help='get report',
                        action='store_true')
    parser.add_argument('-c', '--clean', help='get clean structure',
                        action='store_true')

    parser.add_argument('--get_chain', help='get chain, .e.g A')

    parser.add_argument('--get_seq', help='get seq', action='store_true')

    parser.add_argument('--rosetta2generic', help='convert ROSETTA-like format to a generic pdb',
                        action='store_true')

    parser.add_argument('--get_rnapuzzle_ready', help='get RNApuzzle ready (keep only standard atoms, renumber residues)',
                        action='store_true')

    parser.add_argument('--no_hr', help='do not insert the header into files',
                        action='store_true')

    parser.add_argument('--renumber_residues', help='',
                        action='store_true')

    parser.add_argument('--get_simrna_ready', help='',
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

    args = parser.parse_args()

    s = StrucFile(args.file)
    if args.report:
        print s.get_report()
        print s.get_preview()

    s = StrucFile(args.file)
    if args.clean:
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
        print s.get_text()

    s = StrucFile(args.file)
    if args.get_seq:
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        #print s.get_preview()
        print s.get_seq()

    s = StrucFile(args.file)
    if args.get_chain:
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        #print s.get_preview()
        print s.get_chain(args.get_chain)

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
        print s.get_text()

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
        s.get_simrna_ready(args.renumber_residues)
        print s.get_text()

    if args.get_simrna_ready:
        s = StrucFile(args.file)
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.fix_op_atoms()
        s.renum_atoms()
        if not args.no_hr:
            add_header()
        s.get_simrna_ready(args.renumber_residues)
        print s.get_text()

    if args.renumber_residues:
        s = StrucFile(args.file)
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.get_simrna_ready(args.renumber_residues)
        s.renum_atoms()
        if not args.no_hr:
            add_header()
        print s.get_text()

    if args.delete:
        selection = select_pdb_fragment(args.delete)
        s = StrucFile(args.file)
        if not args.no_hr:
            add_header()
            print 'HEADER --delete ' + args.delete #' '.join(str(selection))
        for l in s.lines:
            if l.startswith('ATOM'):
                chain = l[21]
                resi = int(l[23:26].strip())
                if selection.has_key(chain):
                    if resi in selection[chain]:
                        continue  # print chain, resi
                print l

    if args.edit:
        selection_from, selection_to = select_pdb_fragment(args.edit.split('>')[0]), select_pdb_fragment(args.edit.split('>')[1])
        if len(selection_to) != len(selection_from):
            raise Exception('len(selection_to) != len(selection_from)')
        s = StrucFile(args.file)
        if not args.no_hr:
            add_header()
            print 'HEADER --edit ' + args.edit
        c = 0
        resi_prev = None
        for l in s.lines:
            if l.startswith('ATOM'):
                chain = l[21:22].strip()
                resi = int(l[23:26].strip())
                if selection_from.has_key(chain):
                    if resi in selection_from[chain]:
                        if resi != resi_prev and resi_prev:
                            c += 1
                        resi_prev = resi

                        nl = list(l)
                        nl[21] =  selection_to.keys()[0]
                        nl[23:26] = str(selection_to[selection_to.keys()[0]][c]).rjust(3)
                        nl = ''.join(nl)
                        print nl
                    else:
                        print l

if __name__ == '__main__':
    pass
