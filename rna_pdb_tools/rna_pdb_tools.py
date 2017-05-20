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

    parser.add_argument('--orgmode', help='get a structure in org-mode format <sick!>',
                        action='store_true')

    parser.add_argument('--get_chain', help='get chain, .e.g A')

    parser.add_argument('--fetch', action='store_true', help='fetch file from the PDB db')

    parser.add_argument('--fetch_ba', action='store_true', help='fetch biological assembly from the PDB db')

    parser.add_argument('--get_seq', help='get seq', action='store_true')

    parser.add_argument('--rosetta2generic', help='convert ROSETTA-like format to a generic pdb',
                        action='store_true')

    parser.add_argument('--get_rnapuzzle_ready', help='get RNApuzzle ready (keep only standard atoms, renumber residues)',
                        action='store_true')

    parser.add_argument('--rpr', help='alias to get_rnapuzzle ready)',
                        action='store_true')

    parser.add_argument('--no_hr', help='do not insert the header into files',
                        action='store_true')

    parser.add_argument('--renumber_residues', help='by defult is false',
                        action='store_true')

    parser.add_argument('--collapsed_view', help='',
                        action='store_true')

    parser.add_argument('--cv', help='alias to collapsed_view',
                        action='store_true')

    parser.add_argument('-v', '--verbose', help='tell me more what you\'re doing, please!',
                        action='store_true')

    parser.add_argument('--inplace', help='in place edit the file! [experimental, only for get_rnapuzzle_ready]',
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

# main
if __name__ == '__main__':
    # get version
    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()

    # get parser and arguments
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
            print(add_header(version))
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
            print(add_header(version))
        print(s.get_text())

    if args.get_rnapuzzle_ready or args.rpr:
        if args.inplace:
            shutil.copy(args.file, args.file + '~')
        s = StrucFile(args.file)
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.fix_op_atoms()
        s.renum_atoms()
        s.shift_atom_names()
        s.prune_elements()
        #print s.get_preview()
        #s.write(args.outfile)
        #for l in s.lines:
        #    print l

        remarks = s.get_rnapuzzle_ready(args.renumber_residues, fix_missing_atoms=True, rename_chains=True, verbose=args.verbose)

        if args.inplace:
            with open(args.file, 'w') as f:
                if not args.no_hr:
                    f.write(add_header(version) + '\n')
                if remarks:
                    f.write('\n'.join(remarks) + '\n')
                f.write(s.get_text())

        else:
            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'
            if remarks:    
                output += '\n'.join(remarks) + '\n'
            output += s.get_text()
            try:
                sys.stdout.write(output)
                sys.stdout.flush()
            except IOError:
                pass

    if args.renumber_residues:
        s = StrucFile(args.file)
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.get_rnapuzzle_ready(args.renumber_residues)
        s.renum_atoms()
        if not args.no_hr:
            print(add_header(version))
        print(s.get_text())

    if args.delete:
        if args.inplace:
            shutil.copy(args.file, args.file + '~')

        selection = select_pdb_fragment(args.delete)
        s = StrucFile(args.file)

        output = ''
        if not args.no_hr:
            output += add_header(version) + '\n'
            output += 'HEADER --delete ' + args.delete + '\n' #' '.join(str(selection))
        for l in s.lines:
            if l.startswith('ATOM'):
                chain = l[21]
                resi = int(l[23:26].strip())
                if chain in selection:
                    if resi in selection[chain]:
                        continue  # print chain, resi
                output += l + '\n'

        # write: inplace
        if args.inplace:
            with open(args.file, 'w') as f:
                f.write(output)
        else: # write: to stdout
            try:
                sys.stdout.write(output)
                sys.stdout.flush()
            except IOError:
                pass
                
    if args.edit:
        edit_pdb(args)
        
    if args.fetch:
        fetch(args.file)

    if args.fetch_ba:
        fetch_ba(args.file)

    if args.collapsed_view or args.cv:
        collapsed_view(args)

    if args.orgmode:
        if args.inplace:
            shutil.copy(args.file, args.file + '~')
        s = StrucFile(args.file)
        s.decap_gtp()
        s.fix_resn()
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.fix_op_atoms()
        s.renum_atoms()
        s.shift_atom_names()
        s.prune_elements()
        #print s.get_preview()
        #s.write(args.outfile)
        #for l in s.lines:
        #    print l

        remarks = s.get_rnapuzzle_ready(args.renumber_residues, fix_missing_atoms=True, rename_chains=True, verbose=args.verbose)

        with open(args.file + '~', 'w') as f:
            if not args.no_hr:
                 f.write(add_header(version) + '\n')
            f.write('\n'.join(remarks) + '\n')
            f.write(s.get_text())

        try:
            from Bio import PDB
            from Bio.PDB import PDBIO
            import warnings
            warnings.filterwarnings('ignore', '.*Invalid or missing.*',)
            warnings.filterwarnings('ignore', '.*with given element *',)
        except:
            sys.exit('Error: Install biopython to use this function (pip biopython)')

        parser = PDB.PDBParser()
        struct = parser.get_structure('', args.file + '~')
        model = struct[0]

        # chains [['A', 'seq', [residues]]]
        chains = []
        for c in model.get_list():
            seq = ''
            chain = []
            for r in c:
                chain.append(str(r.get_resname().strip()) + str(r.id[1]))
                seq += r.get_resname().strip()
            chains.append([c.id, seq, chain])

        t = []
        #[['A', 'CCGCCGCGCCAUGCCUGUGGCGG', ['C1', 'C2', 'G3', 'C4', 'C5', 'G6', 'C7', 'G8', 'C9', 'C10', 'A11', 'U12', 'G13', 'C14', 'C15', 'U16', 'G17', 'U18', 'G19', 'G20', 'C21', 'G22', 'G23']], ['B', 'CCGCCGCGCCAUGCCUGUGGCGG', ['C1', 'C2', 'G3', 'C4', 'C5', 'G6', 'C7', 'G8', 'C9', 'C10', 'A11', 'U12', 'G13', 'C14', 'C15', 'U16', 'G17', 'U18', 'G19', 'G20', 'C21', 'G22', 'G23']]]
        for c in chains:
            t.append('* ' + c[0] + ':' + c[2][0][1:] + '-' + c[2][-1][1:] + ' ' + c[1])
            for r in c[2]:
                t.append('** ' + c[0] + ':' + r)
        print '\n'.join(t)

if __name__ == '__main__':
    pass
