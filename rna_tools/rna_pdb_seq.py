#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from icecream import ic
import textwrap
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import os

from rna_tools.rna_tools_lib import edit_pdb, add_header, get_version, \
                          collapsed_view, fetch, fetch_ba, fetch_cif, replace_chain, RNAStructure, \
                          select_pdb_fragment, sort_strings, set_chain_for_struc
from rna_tools.tools.rna_x3dna.rna_x3dna import x3DNA



def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument('--save-to-file', help='<pdb>.fa', action='store_true')

    parser.add_argument('--compact',
                        help=textwrap.dedent("""with --get-seq, get it in compact view'
$ rna_pdb_tools.py --get-seq --compact *.pdb
# 20_Bujnicki_1
ACCCGCAAGGCCGACGGCGCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU # A:1-68
# 20_Bujnicki_2
ACCCGCAAGGCCGACGGCGCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU # A:1-68
# 20_Bujnicki_3
ACCCGCAAGGCCGACGGCGCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU # A:1-68
# 20_Bujnicki_4

"""), action='store_true')

    parser.add_argument('--color-seq', help='color seq, works with --get-seq', action='store_true')

    parser.add_argument('--chain-first', help="", action='store_true')

    parser.add_argument('--fasta',
                        help= textwrap.dedent("""with --get-seq, show sequences in fasta format,
can be combined with --compact (mind, chains will be separated with ' ' in one line)

$ rna_pdb_tools.py --get-seq --fasta --compact input/20_Bujnicki_1.pdb
> 20_Bujnicki_1
ACCCGCAAGGCCGACGGC GCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU

"""), action='store_true')
    parser.add_argument('--oneline', help="", action='store_true')


    parser.add_argument('--uniq', help=textwrap.dedent("""
rna_pdb_tools.py --get-seq --uniq '[:5]' --compact --chain-first * | sort
A:1-121        ACCUUGCGCAACUGGCGAAUCCUGGGGCUGCCGCCGGCAGUACCC...CA # rp13nc3295_min.out.1
A:1-123        ACCUUGCGCGACUGGCGAAUCCUGAAGCUGCUUUGAGCGGCUUCG...AG # rp13cp0016_min.out.1
A:1-123        ACCUUGCGCGACUGGCGAAUCCUGAAGCUGCUUUGAGCGGCUUCG...AG # zcp_6537608a_ALL-000001_AA
A:1-45 57-71   GGGUCGUGACUGGCGAACAGGUGGGAAACCACCGGGGAGCGACCCGCCGCCCGCCUGGGC # solution
"""))

    parser.add_argument('--renum-atoms', help='renumber atoms, tested with --get-seq',
                         action='store_true')
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
        parser = get_parser()
        args = parser.parse_args()

        # quick fix - make a list on the spot
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        analyzed = []
        for f in args.file:
            #####################################
            if args.uniq:
                subname = eval('f' + args.uniq)
                if subname in analyzed:
                    continue
                else:
                    analyzed.append(subname)
            ########
            s = RNAStructure(f)
            if not s.is_pdb():
                print('Error: Not a PDB file %s' % f)
                sys.exit(1)
            s.decap_gtp()
            s.std_resn()
            s.remove_hydrogen()
            s.remove_ion()
            s.remove_water()
            if args.renum_atoms:
                s.renum_atoms()
            s.fix_O_in_UC()
            s.fix_op_atoms()

            output = ''

            # with # is easier to grep this out
            if args.fasta:
                # s.fn vs s.name
                output += s.get_seq(compact=args.compact, chainfirst=args.chain_first, fasta=args.fasta, addfn=s.name, color=args.color_seq) + '\n'
            elif args.oneline:
                output += s.get_seq(compact=args.compact, chainfirst=args.chain_first, color=args.color_seq).strip() + ' # '+ os.path.basename(f.replace('.pdb', '')) + '\n'
            else:
                output += '# ' + os.path.basename(f.replace('.pdb', '')) + '\n'
                output += s.get_seq(compact=args.compact, chainfirst=args.chain_first, color=args.color_seq) + '\n'

            try:
                sys.stdout.write(output)
                sys.stdout.flush()
            except IOError:
                pass

            if args.save_to_file:
                 from rna_pdb_fetch_header import fetch_pdb_header
                 pdb_id = f.split('/')[-1].split('_')[0]
                 header_info = fetch_pdb_header(pdb_id)
                 title = header_info.get('entry', {}).get('struct', {}).get('title', 'Title not found')
                 output = '# ' + title + '\n' + output
                 with open(f.replace('.pdb', '.txt'), 'w') as out:
                    out.write(output)
                    print(output)


