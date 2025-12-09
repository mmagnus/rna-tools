#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_standardize.py - standardzie RNA PDB structures

Usage::

   $ rna_standardize.py <pdb file>

-v is for verbose, --version for version ;)
"""
import argparse
import textwrap
import os
import shutil
import sys
import tempfile
import glob
import os
from rna_tools import rna_cif2pdb
from rna_tools.rna_tools_lib import edit_pdb, add_header, get_version, \
                          collapsed_view, fetch, fetch_ba, replace_chain, RNAStructure, \
                          select_pdb_fragment, sort_strings, set_chain_for_struc
from rna_tools.tools.rna_x3dna.rna_x3dna import x3DNA


def get_parser():
    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', help='', action='version', version=version)

    parser.add_argument('--no-progress-bar', help='for --no-progress-bar for --rpr', action='store_true')

    parser.add_argument('--quiet', help='limit prints', action='store_true')

    parser.add_argument('--renum-nmr', help='',
                        action='store_true')

    parser.add_argument('--inplace', help=textwrap.dedent("""in place edit the file! [experimental,
only for get_rnapuzzle_ready, --delete, --get-ss, --get-seq, --edit-pdb]"""),
                        action='store_true')

    parser.add_argument('-v', '--verbose', help='tell me more what you\'re doing, please!',
                        action='store_true')

    parser.add_argument('--conect', help='add conect to rpr file',  action='store_true')

    parser.add_argument('--conect-no-linkage', help='dont add conect from our residue to another',  action='store_true')

    parser.add_argument('--dont-replace-hetatm', help="replace 'HETATM' with 'ATOM' [tested only with --get-rnapuzzle-ready]",
                        action="store_true")

    parser.add_argument('--keep-hetatm', help='keep hetatoms, [if not replaced anyway with ATOM, see --dont-replace-hetatm', action='store_true')

    parser.add_argument('--here', help=textwrap.dedent("""save a file next to the original file with auto suffix
for --extract it's .extr.pdb"""),
                        action='store_true')

    parser.add_argument('--no-hr', help='do not insert the header into files',
                        action='store_true')

    parser.add_argument('--cif', help='output as cif file',
                        action='store_true')

    parser.add_argument('--check-geometry', help='check connectivity betweeen residues and angles',
                        action='store_true', default=False)

    parser.add_argument('--dont-fix-missing-atoms',
                        help="""used only with --get-rnapuzzle-ready""",
                        action='store_true')

    parser.add_argument('--mdr', help='get structures ready for MD (like rpr but without first)',
                        action='store_true')

    parser.add_argument('--renumber-residues', help='by defult is false',
                        action='store_true')

    parser.add_argument('--suffix', help=textwrap.dedent("""when used with --inplace allows you to change a name of a new file, --suffix del will give <file>_del.pdb (mind added _)"""), default='std')

    parser.add_argument('--dont-report-missing-atoms',
                        help="""used only with --get-rnapuzzle-ready""",
                        action='store_true')

    parser.add_argument('--dont-rename-chains',
                        help=textwrap.dedent("""used only with --get-rnapuzzle-ready.
By default:
   --get-rnapuzzle-ready rename chains from ABC.. to stop behavior switch on this option
"""),
                        action='store_true')

    parser.add_argument('--backbone-only',
                        help="""used only with --get-rnapuzzle-ready, keep only backbone (= remove bases)""",
                        action='store_true')

    parser.add_argument('--no-backbone',
                        help="""used only with --get-rnapuzzle-ready, remove atoms of backbone (define as P OP1 OP2 O5')""",
                        action='store_true')

    parser.add_argument('--bases-only',
                        help="""used only with --get-rnapuzzle-ready, keep only atoms of bases""",
                        action='store_true')

    parser.add_argument('file', help='file', nargs='+')
    #parser.add_argument('outfile', help='outfile')
    return parser


def prepare_input_file(path, version, verbose=False):
    """Return a PDB path for the provided structure, converting CIF when needed."""
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.cif', '.mmcif'):
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='_from_cif.pdb')
        tmp_path = tmp.name
        tmp.close()
        try:
            converted = rna_cif2pdb.convert_cif_to_pdb(
                path,
                add_header_to_output=False,
                version=version,
                verbose=verbose,
                split_conflicting_chains=False,
                output_path=tmp_path,
            )
        except Exception:
            os.remove(tmp_path)
            raise
        if not converted:
            os.remove(tmp_path)
            raise RuntimeError('Failed to convert CIF to PDB: %s' % path)
        if len(converted) > 1:
            os.remove(tmp_path)
            raise RuntimeError('CIF conversion produced multiple PDB files; convert manually and rerun on the desired chain: %s' % path)
        return converted[0], True
    return path, False


# main
if __name__ == '__main__':
    # get version
    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()

    # get parser and arguments
    parser = get_parser()

    args = parser.parse_args()

    remarks_only = False
    # quick fix for one files vs file-s
    if list == type(args.file) and len(args.file) == 1:
        args.file = args.file[0]

    # rnapuzzle
    #import progressbar
    #bar = progressbar.ProgressBar(max_value=len(args.file))
    #bar.update(0)


    if True:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]

        ##################################
        # progress bar only in --inplace mode!
        if args.inplace:
            if not args.no_progress_bar:
                import progressbar
                bar = progressbar.ProgressBar(max_value=len(args.file))
                bar.update(0)

        for c, f in enumerate(args.file):
            if '_std.pdb' in f:
               #print(f'skip {f}')
               continue
            print(f'INPUT: {f} ' + '-' * (80 - len(f)))
            if args.verbose:
                print(f)
            if args.inplace:
                shutil.copy(f, f + '~')

            cleanup_files = []
            try:
                working_file, created_temp = prepare_input_file(f, version, verbose=args.verbose)
                if created_temp:
                    cleanup_files.append(working_file)
                    if args.verbose:
                        print(f' Converted CIF to temporary PDB: {working_file}')

                # keep previous edits
                previous_edits = []
                with open(working_file) as fx:
                    for l in fx:
                        if l.startswith('HEADER --'):
                            previous_edits.append(l.strip())
                ######################
                s = RNAStructure(working_file)
                # ✗ ✓
                if not args.dont_replace_hetatm:
                    if not args.quiet:
                        print('✓ Replace HETATM with ATOM')
                    s.replace_hetatms()

                s.remove_hydrogen()
                if not args.quiet:
                    print('✓ Remove hydrogens')
                s.decap_gtp()
                if not args.quiet:
                    print('✓ Decap GTP')
                s.std_resn()
                s.fix_op_atoms()

                s.remove_ion()
                if not args.quiet:
                    print(f'✓ Remove ions and water')
                s.remove_water()
                # s.renum_atoms()

                s.shift_atom_names()
                s.prune_elements()

                # s.write('tmp.pdb')
                
                rename_chains = False if args.dont_rename_chains else True

                report_missing_atoms = not args.dont_report_missing_atoms
                fix_missing_atom = not args.dont_fix_missing_atoms

                ignore_op3 = False
                if args.mdr:
                    ignore_op3 = True

                remarks = s.get_rnapuzzle_ready(args.renumber_residues, fix_missing_atoms=fix_missing_atom,
                                                rename_chains=rename_chains,
                                                report_missing_atoms=report_missing_atoms,
                                                backbone_only=args.backbone_only,
                                                no_backbone=args.no_backbone,
                                                bases_only=args.bases_only,
                                                keep_hetatm=args.keep_hetatm,
                                                ignore_op3=ignore_op3,
                                                check_geometry=args.check_geometry,
                                            conect=args.conect,
                                            conect_no_linkage=args.conect_no_linkage,
                                                verbose=args.verbose)

                output = ''
                if not args.no_hr:
                    output += add_header(version) + '\n'
                if remarks:
                    output += '\n'.join(remarks) + '\n'
                output += s.get_text() + '\n'

                if remarks_only:
                    sys.stdout.write('\n'.join(remarks))
                    sys.stdout.flush()
                else:
                    args.here = True
                    if args.here:
                        fileroot, _ = os.path.splitext(f)
                        nf = f'{fileroot}_{args.suffix}.pdb'
                        with open(nf, 'w') as fio:
                                for r in remarks:
                                    print(r.replace('REMARK 250', ' '))
                                if not args.quiet: print('Output:', nf)
                                print(s.get_seq(addfn = nf))                                
                                fio.write(output)

                        if args.cif:
                            cif = f'{fileroot}_{args.suffix}.cif'
                            from Bio.PDB import PDBParser, MMCIFIO
                            pdb_parser = PDBParser()
                            structure = pdb_parser.get_structure('structure_id', nf)
                            mmcif_writer = MMCIFIO()
                            mmcif_writer.set_structure(structure)
                            mmcif_writer.save(cif)
                            print(f'Output: {cif}')
                            os.remove(nf)
                    else:
                        try:
                            sys.stdout.write(output)
                            sys.stdout.flush()
                        except IOError:
                            pass
            finally:
                for tmp in cleanup_files:
                    if os.path.exists(tmp):
                        os.remove(tmp)
            # bar.update(c)

        # hmm... fix for problem with renumbering, i do renumbering
        # and i stop here
        # i'm not sure that this is perfect
