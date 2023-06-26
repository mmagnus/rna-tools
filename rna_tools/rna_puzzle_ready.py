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

    parser.add_argument('--keep-hetatm', help='keep hetatoms', action='store_true')

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

    parser.add_argument('--replace-hetatm', help="replace 'HETATM' with 'ATOM' [tested only with --get-rnapuzzle-ready]",
                        action="store_true")

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
            if args.replace_hetatm:
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
                # progress bar only in --inplace mode!
                if not args.no_progress_bar:
                   bar.update(c)

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
                        if '_rpr' not in f:  # good idea?
                            nf = f.replace('.pdb', '_rpr.pdb')
                            with open(nf, 'w') as fio:
                                print(nf)
                                fio.write(output)
                    else:
                        try:
                            sys.stdout.write(output)
                            sys.stdout.flush()
                        except IOError:
                            pass
        # hmm... fix for problem with renumbering, i do renumbering
        # and i stop here
        # i'm not sure that this is perfect
        sys.exit(0)


    if args.renumber_residues:
        s = RNAStructure(args.file)
        s.remove_hydrogen()
        s.remove_ion()
        s.remove_water()
        s.get_rnapuzzle_ready(args.renumber_residues)
        s.renum_atoms()
        if not args.no_hr:
            print(add_header(version))
        print(s.get_text())

# args.renumber_resides_dirty
    if args.renum_residues_dirty:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')

            s = RNAStructure(f)

            output = ''
            #if not args.no_hr:
            #    output += add_header(version) + '\n'
            #    output += 'HEADER --delete ' + args.delete + '\n'  # ' '.join(str(selection))
            c = 1
            old_resi = -1
            for l in s.lines:
                if l.startswith('ATOM') or l.startswith("HETATOM"):
                    resi = int(l[23:26].strip())
                    if resi != old_resi:
                        old_resi = resi
                        c += 1
                    # print(resi, c)
                    #resi = c
                    #if chain in selection:
                    #    if resi in selection[chain]:
                    #        continue  # print chain, resi
                    output += l[:23] + str(c).rjust(3) + l[26:] + '\n'
            # write: inplace
            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass


    if args.undo:
        # quick fix - make a list on the spot
        dir = args.file #os.path.abso(os.path.dirname(args.file))
        for f in glob.glob(dir + '/*.pdb~'):
            if args.verbose:
                print(f, '->',f.replace('.pdb~', '.pdb'))
            os.rename(f, f.replace('.pdb~', '.pdb'))

    if args.delete:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')

            selection = select_pdb_fragment(args.delete)
            s = RNAStructure(f)

            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'
                output += 'HEADER --delete ' + args.delete + '\n'  # ' '.join(str(selection))
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
                os.rename(f, f.replace('.pdb', '.pdb~'))                
                if args.suffix:
                    f = f.replace('.pdb', '_' + args.suffix + '.pdb')
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass

    if args.replace_chain:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')

            # --replace_chain <file> without A:<file> it will be easier than --x "A:<file>"
            s = RNAStructure(args.replace_chain)
            chain_ids = (s.get_all_chain_ids())
            if len(chain_ids) > 1:
                raise Exception('There is more than one chain in the inserted PDB file. There should be only one chain, the one you want to insert to the PDB.')
            out = replace_chain(f, args.replace_chain, list(chain_ids)[0])
            print(out)

    if args.mutate:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        from rna_tools.tools.mini_moderna3.moderna import *

        for f in args.file:
            if args.ignore_files:
                if args.ignore_files in f:
                    continue
                
            if args.inplace:
                shutil.copy(f, f + '~')  # create a backup copy if inplace

            # create a working copy of the main file
            ftf = tempfile.NamedTemporaryFile(delete=False).name  # f temp file
            shutil.copy(f, ftf)  # create a backup copy if inplace

            # go over each chain
            # rna_pdb_tools.py --mutate 'A:1CB:1G,A:1U+B:1A' CG_AB.pdb > ~/Desktop/a.pdb
            args.mutate = args.mutate.upper()
            for m in args.mutate.split(','):  # A:1A, B:1A
                chain, resi_mutate_to = m.strip().split(':')  # A:1A+2C
                resi_mutate_to_list = resi_mutate_to.split('+')  # A:1A+2C

                model = load_model(f, chain)
                # go over each mutation in this chain
                for resi_mutate_to in resi_mutate_to_list:
                    resi = resi_mutate_to[:-1]
                    mutate_to = resi_mutate_to[-1]
                    if args.verbose:
                        print('Mutate', model[resi], 'to', mutate_to)
                    exchange_single_base(model[resi], mutate_to)

                # multi mutation in one chain
                tf = tempfile.NamedTemporaryFile(delete=False)
                model.write_pdb_file(tf.name)

                # work on the copy of f, ftf
                output = replace_chain(ftf, tf.name, chain)
                with open(ftf, 'w') as tmp:
                    tmp.write(output)

            # write: inplace
            if args.inplace:
                # ftf now is f, get ready for the final output
                if args.suffix:
                    f = f.replace('.pdb', '_' + args.suffix + '.pdb')
                # rpr on the file
                shutil.copy(ftf, f)
                os.system('rna_pdb_tools.py --rpr --no-progress-bar --inplace ' + f)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass


# extract
    if args.extract:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')

            selection = select_pdb_fragment(args.extract)
            s = RNAStructure(f)

            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'
                output += 'HEADER extract ' + args.extract + '\n'  # ' '.join(str(selection))
            for l in s.lines:
                if l.startswith('ATOM'):
                    chain = l[21]
                    resi = int(l[23:26].strip())
                    if chain in selection:
                        if resi in selection[chain]:
                            # continue  # print chain, resi
                            output += l + '\n'

            # write: inplace
            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            elif args.here:
                if '_extr' not in f:  # good idea?
                    nf = f.replace('.pdb', '_extr.pdb')
                    with open(nf, 'w') as fio:
                        print(nf)
                        fio.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass


    if args.extract_chain:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')

            selection = select_pdb_fragment(args.extract)
            s = RNAStructure(f)

            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'
                output += 'HEADER extract ' + args.extract + '\n'  # ' '.join(str(selection))
            for l in s.lines:
                if l.startswith('ATOM') or l.startswith('TER') or l.startswith('HETATM'):
                    chain = l[21]
                    if chain in args.extract_chain:
                        output += l + '\n'

            # write: inplace
            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass

    if args.is_pdb:
        s = RNAStructure(args.file)
        output = str(s.is_pdb())
        sys.stdout.write(output + '\n')

    if args.un_nmr:
        s = RNAStructure(args.file)
        str(s.un_nmr())

    if args.is_nmr:
        struc = RNAStructure(args.file)
        output = str(struc.is_nmr(args.verbose))
        sys.stdout.write(output + '\n')

    #edit
    if args.edit:
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.verbose:
                print(f)
            if args.inplace:
                shutil.copy(f, f + '~')

            output = edit_pdb(f, args)

            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass


    if args.fetch:
        fn = fetch(args.file)
            
    if args.fetch_chain:
        fn = fetch(args.file[0])
        chain = args.file[1]
        nfn = fn.replace('.pdb','') + '_' + chain + '.pdb' # 6bk8_H.pdb
        #cmd = 'rna_pdb_tools.py --get-chain %s %s' > %s' % (chain, fn, nfn)
        cmd = 'rna_pdb_tools.py --get-chain %s %s' % (chain, fn)
        os.system(cmd)
        # print(nfn)

    if args.fetch_ba:
        fetch_ba(args.file)

    if args.collapsed_view or args.cv:
        collapsed_view(args)

    if args.delete_anisou:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')

            s = RNAStructure(f)

            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'

            for l in s.lines:
                if l.startswith('ANISOU'):
                    continue
                else:
                    output += l + '\n'

            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass


    if args.swap_chains:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')
            # rename_chain 'A>B'
            s = RNAStructure(f)
            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'

            chain_id_old, chain_id_new = args.swap_chains.split('>')
            if chain_id_old == chain_id_new:
                output = open(f).read() # return the file ;-) itself unchanged
            else:
                s.rename_chain(chain_id_new, '_')
                s.rename_chain(chain_id_old, chain_id_new)
                output += s.rename_chain('_', chain_id_old)

            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass

    if args.rename_chain:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')
            # rename_chain 'A>B'
            s = RNAStructure(f)
            chain_id_old, chain_id_new = args.rename_chain.split('>')
            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'
            output += s.rename_chain(chain_id_old, chain_id_new)
            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass


    if args.split_alt_locations:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            if args.inplace:
                shutil.copy(f, f + '~')

            #s = RNAStructure(f)
            s = open(f)
            output = ''
            if not args.no_hr:
                output += add_header(version) + '\n'

            # First, collect all alternative locations.
            alts = set()
            for l in s:
                if l.startswith('ATOM'):
                    if l[16].strip():
                        alts.add(l[16])
            s.close()

            if args.verbose:
                print('alts:', alts)

            for index, alt in enumerate(alts):
                output += 'MODEL %i' % (index + 1)
                s = open(f)
                for l in s:
                    if l.startswith('ATOM') or l.startswith('HETATM'):
                        # if empty, then print this line
                        if l[16] == ' ':
                            output += l
                        if l[16] == alt:
                            output += l
                    else:
                        output += l
                output += 'ENDMDL\n'
                s.close()

            if args.inplace:
                with open(f, 'w') as f:
                    f.write(output)
            else:  # write: to stdout
                try:
                    sys.stdout.write(output)
                    sys.stdout.flush()
                except IOError:
                    pass

    if args.fix:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            cmd = 'pdbfixer ' + f + ' --add-atoms all --add-residues'
            print(cmd)
            os.system(cmd)
            if args.inplace:
                shutil.move("output.pdb", f)
            else:
                shutil.move("output.pdb", f.replace('.pdb', '_fx.pdb'))
                            
    if args.to_mol2:
        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            cmd = 'obabel -i pdb ' + f + ' -o mol2 -O ' + f.replace('.pdb', '.mol2')
            print(cmd)
            os.system(cmd)


    if args.nmr_dir:
        files = sort_strings(glob.glob(args.nmr_dir + '/' + args.file))

        c = 1

        for f in files:
            #s = RNAStructure(f)
            print("MODEL        " + str(c))
            # at some point I could use RNAStructure for this
            print(open(f).read())
            #print(s.get_text(add_end=False))
            print('ENDMDL')
            c += 1
        print('END')

    if args.set_chain:
        if list != type(args.file):
            args.file = [args.file]

        for f in args.file:
            txt = set_chain_for_struc(f, args.set_chain)
            if args.inplace:
                shutil.move(f, f.replace('.pdb', '.pdb~'))
                with open(f, 'w') as fi:
                    fi.write(txt)
            else:
                print(txt)

    if args.replace_htm:
        if list != type(args.file):
            args.file = [args.file]

        for f in args.file:
            with open(f) as fn:
                txt = fn.read()
            txt = txt.replace('HETATM', 'ATOM  ')
            if args.inplace:
                shutil.move(f, f.replace('.pdb', '.pdb~'))
                with open(f, 'w') as fi:
                    fi.write(txt)
            else:
                print(txt)

    if args.renum_nmr:
        if list != type(args.file):
            args.file = [args.file]
        txt = ''
        for f in args.file:
            c = 1
            for l in open(f):
                if l.startswith('MODEL'):
                    txt += "MODEL       " + str(c) + '\n'
                    c += 1
                elif l.strip() == 'END':
                    pass
                else:
                    txt += l
        if args.inplace:
            shutil.move(f, f.replace('.pdb', '.pdb~'))
            with open(f) as fi:
                fi.write(txt)
        else:
            print(txt)
            
    from rna_tools.rna_tools_config import PYMOL_PATH
    sys.path.insert(0, PYMOL_PATH)
    if args.cif2pdb:
        try:
            from pymol import cmd
        except ImportError:
            print('This functionality needs PyMOL. Install it or if installed, check your setup')
            sys.exit(1)

        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            cmd.load(f)
            fo = f.replace('.cif', '.pdb')
            cmd.save(fo, '(all)')
            cmd.delete('all')
            print(fo, 'saved')
            
    if args.pdb2cif:
        try:
            from pymol import cmd
        except ImportError:
            print('This functionality needs PyMOL. Install it or if installed, check your setup')
            sys.exit(1)

        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for f in args.file:
            cmd.load(f)
            fo = f.replace('.pdb', '.cif')
            cmd.save(fo, '(all)')
            cmd.delete('all')
            print(fo, 'saved')
