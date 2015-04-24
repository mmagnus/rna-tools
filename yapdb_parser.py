#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import os
import time

from pdb_parser_lib import *

def add_header():
    now = time.strftime("%c")
    print 'HEADER Generated with yapdb_parser, version: %s (https://github.com/m4rx9/rna-pdb-tools) %s' % (version, now)

if __name__ == '__main__':
    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()
    parser = argparse.ArgumentParser('yapdb_parser ver: %s' % version)
    parser.add_argument('-r', '--report', help='get report',
                        action='store_true')
    parser.add_argument('-c', '--clean', help='get clean structure',
                        action='store_true')

    parser.add_argument('--getchain', help='get chain, .e.g A')

    parser.add_argument('--getseq', help='get seq', action='store_true')

    parser.add_argument('--rosetta2generic', help='convert ROSETTA-like format to generic pdb',
                        action='store_true')

    parser.add_argument('--getrnapuzzle', help='get RNApuzzle ready',
                        action='store_true')

    parser.add_argument('--nohr', help='do not insert the header into files',
                        action='store_true')

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
        if not args.nohr:
            add_header()
        print s.get_text()

    s = StrucFile(args.file)
    if args.getseq:
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
        #s.write(args.outfile)

    s = StrucFile(args.file)
    if args.getchain:
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.renum_atoms()
        s.fix_O_in_UC()
        s.fix_op_atoms()
        #print s.get_preview()
        print s.get_chain(args.getchain)
        #s.write(args.outfile)


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
        if not args.nohr:
            add_header()
        print s.get_text()

    if args.getrnapuzzle:
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
        if not args.nohr:
            add_header()
        s.get_rnapuzzle_ready()
        print s.get_text()
