#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""rna_pdb_tools - a swiss army knife to manipulation of RNA pdb structures

Usage::

   $ rna_pdb_tools.py *

See `rna_pdb_merge_into_one.py` to merge PDB files in the order as you like into one NMR-like (multimodel) file

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

    parser.add_argument('--renum-nmr', help='',
                        action='store_true')

    parser.add_argument('--inplace', help=textwrap.dedent("""in place edit the file! [experimental,
only for get_rnapuzzle_ready, --delete, --get-ss, --get-seq, --edit-pdb]"""),
                        action='store_true')

    parser.add_argument('-v', '--verbose', help='tell me more what you\'re doing, please!',
                        action='store_true')

    parser.add_argument('--dont-replace-hetatm', help="replace 'HETATM' with 'ATOM' [tested only with --get-rnapuzzle-ready]",
                        action="store_true")

    parser.add_argument('--keep-hetatm', help='keep hetatoms, [if not replaced anyway with ATOM, see --dont-replace-hetatm', action='store_true')

    parser.add_argument('--here', help=textwrap.dedent("""save a file next to the original file with auto suffix
for --extract it's .extr.pdb"""),
                        action='store_true')

    parser.add_argument('--no-hr', help='do not insert the header into files',
                        action='store_true')

    parser.add_argument('--dont-fix-missing-atoms',
                        help="""used only with --get-rnapuzzle-ready""",
                        action='store_true')

    parser.add_argument('--mdr', help='get structures ready for MD (like rpr but without first)',
                        action='store_true')

    parser.add_argument('--renumber-residues', help='by defult is false',
                        action='store_true')

    parser.add_argument('--suffix', help=textwrap.dedent("""when used with --inplace allows you to change a name of a new file, --suffix del will give <file>_del.pdb (mind added _)"""))

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
            if args.verbose:
                print(f)
            if args.inplace:
                shutil.copy(f, f + '~')

            # keep previous edits
            previous_edits = []
            with open(f) as fx:
                for l in fx:
                    if l.startswith('HEADER --'):
                        previous_edits.append(l.strip())
            ######################

            s = RNAStructure(f)
            if not args.dont_replace_hetatm:
                s.replace_hetatms()

            s.remove_hydrogen()
            s.decap_gtp()
            s.std_resn()
            s.fix_op_atoms()

            s.remove_ion()
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
                                            verbose=args.verbose)


            if args.inplace:
                if args.suffix:
                    f = f.replace('.pdb', '_' + args.suffix + '.pdb')
                with open(f, 'w') as f:
                    if not args.no_hr:
                        f.write(add_header(version) + '\n')
                    if previous_edits:
                        f.write('\n'.join(previous_edits) + '\n')
                    if remarks:
                        f.write('\n'.join(remarks) + '\n')
                    f.write(s.get_text())

            else:
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
                            nf = f.replace('.pdb', '_std.pdb')
                            with open(nf, 'w') as fio:
                                print('Output:', nf)
                                print(s.get_seq())                                
                                fio.write(output)


                    else:
                        try:
                            sys.stdout.write(output)
                            sys.stdout.flush()
                        except IOError:
                            pass
            # bar.update(c)

        # hmm... fix for problem with renumbering, i do renumbering
        # and i stop here
        # i'm not sure that this is perfect
