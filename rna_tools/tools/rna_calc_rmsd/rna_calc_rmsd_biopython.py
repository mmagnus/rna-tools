#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Examples::

  rna_calc_rmsd_biopython.py -t 2nd_triplex_FB_ABC.pdb  triples-all-v2-rpr/*.pdb --way backbone+sugar --save --suffix ABC --result abc.csv

  rna_calc_rmsd_biopython.py -t 2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples2/Triple_cWW_tSH_GCA_exemplar_rpr_ren.pdb  --way backbone+sugar --save --column-name 'rmsd' --triple-mode > out.txt

"""
from __future__ import print_function

import Bio.PDB.PDBParser
import Bio.PDB.Superimposer
from Bio.Align import PairwiseAligner
import copy
import warnings
warnings.filterwarnings('ignore', '.*Invalid or missing.*',)
warnings.filterwarnings('ignore', '.*with given element *',)
warnings.filterwarnings('ignore', '.*is discontinuous*',)
warnings.filterwarnings('ignore', '.*Ignoring unrecognized record*',)

from collections import OrderedDict
from rna_tools.rna_tools_lib import set_chain_for_struc, RNAStructure
import argparse
import sys
import glob
import re
import os

NUCLEOTIDE_THREE_TO_ONE = {
    'A': 'A', 'DA': 'A', 'ADE': 'A',
    'C': 'C', 'DC': 'C', 'CYT': 'C',
    'G': 'G', 'DG': 'G', 'GUA': 'G',
    'U': 'U', 'DU': 'U', 'DT': 'U', 'T': 'U', 'PSU': 'U', 'URA': 'U',
    'I': 'I', 'DI': 'I'
}

WAY_TO_ATOMS = {
    'c1': ["C1'"],
    'backbone+sugar': "P,OP1,OP2,C5',O5',C4',O4',C3',O3',C2',O2',C1'".split(',')
}


def _residue_is_polymer(residue):
    hetfield = residue.id[0].strip()
    return hetfield == ''


def _residue_to_letter(residue):
    name = residue.get_resname().strip().upper()
    return NUCLEOTIDE_THREE_TO_ONE.get(name, 'N')


def _build_alignment_strings(alignment, seq1, seq2):
    seq1_line = []
    seq2_line = []
    match_line = []
    matched_pairs = []
    matches = 0

    path = getattr(alignment, 'path', None)
    if not path:
        return '', '', '', matched_pairs, matches

    seq1_idx = path[0][0]
    seq2_idx = path[0][1]

    for next_i, next_j in path[1:]:
        while seq1_idx < next_i or seq2_idx < next_j:
            if seq1_idx < next_i and seq2_idx < next_j:
                char_a = seq1[seq1_idx]
                char_b = seq2[seq2_idx]
                seq1_line.append(char_a)
                seq2_line.append(char_b)
                if char_a == char_b:
                    match_line.append('|')
                    matches += 1
                else:
                    match_line.append(' ')
                matched_pairs.append((seq1_idx, seq2_idx))
                seq1_idx += 1
                seq2_idx += 1
            elif seq1_idx < next_i:
                char_a = seq1[seq1_idx]
                seq1_line.append(char_a)
                seq2_line.append('-')
                match_line.append(' ')
                seq1_idx += 1
            else:
                char_b = seq2[seq2_idx]
                seq1_line.append('-')
                seq2_line.append(char_b)
                match_line.append(' ')
                seq2_idx += 1

    seq1_line = ''.join(seq1_line)
    match_line = ''.join(match_line)
    seq2_line = ''.join(seq2_line)
    return seq1_line, match_line, seq2_line, matched_pairs, matches


def _print_alignment(seq1_line, match_line, seq2_line, header,
                     target_len=None, model_len=None, matches=None, residue_pairs=None):
    print(header)
    if seq1_line:
        print(seq1_line)
    if match_line:
        print(match_line)
    if seq2_line:
        print(seq2_line)
    details = []
    if target_len is not None:
        details.append('target_len={}'.format(target_len))
    if model_len is not None:
        details.append('model_len={}'.format(model_len))
    if residue_pairs is not None:
        details.append('residue_pairs={}'.format(residue_pairs))
    if matches is not None:
        details.append('matches={}'.format(matches))
    if details:
        print('# alignment info:', ', '.join(details))


class RNAmodel:
    def __init__(self, fpath):
        # parser 1-5 -> 1 2 3 4 5
        self.struc = Bio.PDB.PDBParser().get_structure('', fpath)
        self.__get_atoms()
        self.fpath = fpath
        self.fn = os.path.basename(fpath)
        # self.atoms = []
        # if save:
        #    self.save() # @save

    def __get_atoms(self):
        self.atoms = []
        self.residues = []
        self.sequence = ''
        self.sequence_residues = []
        self.sequence_positions = []
        self.sequence_atom_maps = []
        for res in self.struc.get_residues():
            self.residues.append(res)
            atom_map = OrderedDict()
            for at in res:
                self.atoms.append(at)
                atom_map[at.name] = at
            if _residue_is_polymer(res):
                self.sequence += _residue_to_letter(res)
                chain_id = res.get_parent().id
                resseq = res.get_id()[1]
                self.sequence_positions.append(f"{chain_id}:{resseq}")
                self.sequence_residues.append(res)
                self.sequence_atom_maps.append(atom_map)
        return self.atoms

    def __str__(self):
        return self.fn #+ ' # beads' + str(len(self.residues))

    def __repr__(self):
        return self.fn #+ ' # beads' + str(len(self.residues))

    def get_report(self):
        """Str a short report about rna model""" 
        t = ' '.join(['File: ', self.fn, ' # of atoms:', str(len(self.atoms)), '\n'])
        for r,a in zip(self.residues, self.atoms ):
            t += ' '.join(['resi: ', str(r) ,' atom: ', str(a) , '\n' ])
        return t

    def _atoms_from_sequence_alignment(self, other_rnamodel, way):
        if not self.sequence or not other_rnamodel.sequence:
            raise ValueError('Cannot align sequences when one of the models has no polymer residues')

        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        gap_open_penalty = -abs(args.align_gap_open if args.align_gap_open is not None else 0.0)
        gap_extend_source = args.align_gap_extend if args.align_gap_extend is not None else args.align_gap_open
        gap_extend_penalty = -abs(gap_extend_source if gap_extend_source is not None else 0.0)
        aligner.open_gap_score = gap_open_penalty
        aligner.extend_gap_score = gap_extend_penalty

        alignments = aligner.align(self.sequence, other_rnamodel.sequence)
        try:
            alignment = alignments[0]
        except IndexError:
            raise ValueError('Biopython did not return an alignment')

        seq1_line, match_line, seq2_line, matched_pairs, matches = _build_alignment_strings(
            alignment, self.sequence, other_rnamodel.sequence)

        if getattr(args, 'print_alignment', False):
            header = '# alignment between {} and {}'.format(self.fn, other_rnamodel.fn)
            _print_alignment(
                seq1_line,
                match_line,
                seq2_line,
                header,
                target_len=len(seq1_line),
                model_len=len(seq2_line),
                matches=matches,
                residue_pairs=len(matched_pairs))

        if not matched_pairs:
            raise ValueError('Sequence alignment returned no overlapping residues')

        allowed_names = WAY_TO_ATOMS.get(way)
        atoms_self = []
        atoms_other = []
        residues_used = 0

        for idx_a, idx_b in matched_pairs:
            atom_map_a = self.sequence_atom_maps[idx_a]
            atom_map_b = other_rnamodel.sequence_atom_maps[idx_b]
            if allowed_names:
                atom_names = [n for n in allowed_names if n in atom_map_a and n in atom_map_b]
            else:
                atom_names = [n for n in atom_map_a.keys() if n in atom_map_b]
            if not atom_names:
                continue
            for name in atom_names:
                atoms_self.append(atom_map_a[name])
                atoms_other.append(atom_map_b[name])
            residues_used += 1

        if not atoms_self:
            raise ValueError('No overlapping atoms remained after sequence alignment')

        seq_len = max(len(self.sequence), len(other_rnamodel.sequence))
        identity = (matches / float(seq_len)) if seq_len else 0.0

        info = {
            'aligned_residues': residues_used,
            'identity': identity,
            'atoms': len(atoms_self),
            'alignment_pairs': len(matched_pairs)
        }
        return atoms_self, atoms_other, info

    def get_rmsd_to(self, other_rnamodel, way="", triple_mode=False, save=True):
        """Calc rmsd P-atom based rmsd to other rna model"""
        sup = Bio.PDB.Superimposer()
        if args.align_sequence:
            if triple_mode:
                raise ValueError('Sequence alignment mode is not supported together with triple mode')
            self.atoms_for_rmsd, other_atoms_for_rmsd, info = self._atoms_from_sequence_alignment(other_rnamodel, way)
            if args.debug:
                print('Aligned residues:', info['aligned_residues'], 'identity:', round(info['identity'], 3), '#atoms:', info['atoms'])
            if not self.atoms_for_rmsd:
                raise ValueError('No overlapping atoms remained after sequence alignment based filtering')
        elif way == 'c1':
            self.atoms_for_rmsd = []
            for a in self.atoms:
                if a.name in ["C1'"]:
                    self.atoms_for_rmsd.append(a)
            if args.debug: print('atoms_for_rmsd', len(self.atoms_for_rmsd))

            other_atoms_for_rmsd = []
            for a in other_rnamodel.atoms:
                if a.name in ["C1'"]:
                    other_atoms_for_rmsd.append(a)
            if args.debug: print('other_atoms_for_rmsd', len(other_atoms_for_rmsd))

        elif way == 'backbone+sugar':
            self.atoms_for_rmsd = []
            for a in self.atoms:
                if a.name in "P,OP1,OP2,C5',O5',C4',O4',C3',O3',C2',O2',C1'".split(','):
                    self.atoms_for_rmsd.append(a)
            if args.debug: print('atoms_for_rmsd', len(self.atoms_for_rmsd))
                                 
            other_atoms_for_rmsd = []

            for a in other_rnamodel.atoms:
                if a.name in "P,OP1,OP2,C5',O5',C4',O4',C3',O3',C2',O2',C1'".split(','):
                    other_atoms_for_rmsd.append(a)
            if args.debug: print('other_atoms_for_rmsd', len(other_atoms_for_rmsd))
        else:
            self.atoms_for_rmsd = self.atoms
            other_atoms_for_rmsd = other_rnamodel.atoms

        if triple_mode:
            def chunks(lst, n):
                """Yield successive n-sized chunks from lst.
                https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
                """
                for i in range(0, len(lst), n):
                    yield lst[i:i + n]

            rmsd_min = 10000 # ugly
            import itertools
            per = list(itertools.permutations([0, 1, 2]))
            lst = list(chunks(other_atoms_for_rmsd,
                              int(len(other_atoms_for_rmsd)/3))) # for 1 atom, this will be 1 x 3 [3 residues]
                       # so so len is 3 atoms so / by 3 to get how many atoms per residue
            sup_min = None

            for p in per:
                patoms = []
                for i in p: # p=(1, 2, 3)
                    patoms.extend(lst[i])
                #print(self.atoms_for_rmsd)
                ## print('patoms')
                ## for a in patoms:
                ##     print(a, a.get_parent().get_id())
                ## print('self.atoms_for_rmsd')
                ## for a in self.atoms_for_rmsd:
                ##     print(a, a.get_parent().get_id())
                #sup.set_atoms(patoms, self.atoms_for_rmsd)
                sup.set_atoms(self.atoms_for_rmsd, patoms)
                rms = round(sup.rms, 3)

                rt = None
                seq = ''
                for a in patoms:
                    r = a.get_parent()
                    if r != rt:
                        rt = r
                        seq += r.get_resname().strip()

                if rms < rmsd_min:
                    rmsd_min = rms
                    sup_min = copy.copy(sup)
                    suffix = seq
                    p_min = p
                    seq_min = seq

                if args.debug:
                    print(p, '', [i + 1 for i in p], end=' ')
                    print(seq, 'seq_min: a' + seq_min, rms)

            index = [0 ,0 ,0]
            index[0] = p_min.index(0)
            index[1] = p_min.index(1)
            index[2] = p_min.index(2)

            # ugly re-set 123 to crazy id ! + 100, so this will
            # fill up 1 2 3 for the second for
            rs = other_rnamodel.struc[0].get_residues()
            for i, r in enumerate(rs):
                r.id = (' ', index[i] + 153, ' ') # ugly, some random offset
            
            for i, r in enumerate(other_rnamodel.struc[0].get_residues()):
                r.id = (' ', index[i] + 1, ' ')
                if args.debug: print(r)

            io = Bio.PDB.PDBIO()
            sup_min.apply(other_rnamodel.struc.get_atoms())
            # if args.debug: print(p_min, [i + 1 for i in p_min])
            io.set_structure(other_rnamodel.struc)
            
            args.save_here = True
            if args.save_here:
                f = os.path.basename(self.fpath.replace('.pdb', '_aligned'))
                # print(f)
                try:
                    os.mkdir(f)
                except:
                    pass
                fout = f + os.sep + os.path.basename(other_rnamodel.fpath.replace('.pdb', '_aligned_s' + suffix + '.pdb'))
            else:
                fout = other_rnamodel.fpath.replace('.pdb', '_aligned_s' + suffix + '.pdb')

            # if args.debug: print(fout)
            io.save(fout)

            # ugly set chain to A
            set_chain_for_struc(fout, 'A', save_file_inplace=True)
            # and now run this to sort into 1 2 3

            r = RNAStructure(fout)
            remarks = r.get_remarks_text()
            r1 = r.get_res_text('A', 1)
            r2 = r.get_res_text('A', 2)
            r3 = r.get_res_text('A', 3)
            with open(fout, 'w') as f:
                f.write(remarks)
                f.write(r1)
                f.write(r2)
                f.write(r3)
            r.reload()
            r.get_rnapuzzle_ready()
            r.write()
            return str(rmsd_min) + ',s' + seq_min + ',' + os.path.basename(fout)

        else:
            sup.set_atoms(self.atoms_for_rmsd, other_atoms_for_rmsd)
            rms = round(sup.rms, 2)

            if save:
                ## io = Bio.PDB.PDBIO()
                ## sup.apply(self.struc.get_atoms())
                ## io.set_structure(self.struc)
                ## io.save("aligned.pdb")

                io = Bio.PDB.PDBIO()
                sup.apply(other_rnamodel.struc.get_atoms())
                io.set_structure(other_rnamodel.struc)
                io.save(other_rnamodel.fpath.replace('.pdb', '_' + args.suffix + '.pdb'))
        return rms

    
def get_rna_models_from_dir(directory):
    models = []
    if not os.path.exists(directory):
        raise Exception('Dir does not exist! ', directory)
    files = glob.glob(directory + "/*.pdb")
    files_sorted = sort_nicely(files)
    for f in files_sorted:
        #print f
        models.append(RNAmodel(f))
    return models

def sort_nicely( l ):
   """ Sort the given list in the way that humans expect.

   http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
   """
   convert = lambda text: int(text) if text.isdigit() else text
   alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
   l.sort( key=alphanum_key )
   return l


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t',"--target",
                         help="target file")
    parser.add_argument('--result', #default='out.rmsd',
                         help="result file")
    parser.add_argument('--ignore-files', default='aligned', help="use to ignore files, .e.g. with 'aligned'")
    parser.add_argument('--suffix', default='aligned', help="used with --saved, by default: aligned")
    parser.add_argument('--way', help="e.g., backbone+sugar")
    parser.add_argument('--align-sequence', action='store_true',
                        help='align sequences with Biopython and compute RMSD only over aligned residues')
    parser.add_argument('--align-gap-open', type=float, default=10.0,
                        help='gap opening penalty (positive value) when using --align-sequence, default: 10')
    parser.add_argument('--align-gap-extend', type=float,
                        help='gap extension penalty (positive value) when using --align-sequence; default: same as gap open')
    parser.add_argument('--print-alignment', action='store_true',
                        help='print the pairwise sequence alignment used for atom matching')
    parser.add_argument('--print-target-sequence', action='store_true',
                        help='print the polymer sequence extracted from the target structure and continue processing')
    parser.add_argument('--triple-mode', help="same crazy strategy to align triples", action="store_true")
    parser.add_argument('--column-name', help="name column for rmsd, by default 'rmsd', but can also be a name of the target file",
                        default="rmsd")
    parser.add_argument("-s", "--save", action="store_true", help="set suffix with --suffix, by default: aligned")
    parser.add_argument('files', help='files', nargs='+')
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--sort", action="store_true", help='sort results based on rmsd (ascending)')
    return parser

# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    target = RNAmodel(args.target)

    if args.print_target_sequence:
        print('# target sequence (polymer residues only):', target.fn)
        if target.sequence:
            print(target.sequence)
            print('# length:', len(target.sequence))
        else:
            print('# warning: target structure does not contain polymer residues with standard names')

    models = args.files  # opts.input_dir
    tmp = []
    if args.ignore_files:
        for f in args.files:
            if args.debug: print(f)
            if args.ignore_files in f:
                continue
            tmp.append(f)
        models = tmp
 
    print('# of models:', len(models))
    c = 1
    t = 'fn,' + args.column_name + ',aligned_seq, aligned_fn\n'
    for m in models:
        mrna = RNAmodel(m)
        #print r1.fn, r2.fn, r1.get_rmsd_to(r2)#, 'tmp.pdb')
        # print(m)
        try:
            rmsd = target.get_rmsd_to(mrna, args.way, args.triple_mode, args.save) #, 'tmp.pdb')
        except ValueError as exc:
            if args.align_sequence:
                print('# warning:', mrna.fn, exc, file=sys.stderr)
                rmsd = 'nan'
            else:
                raise
        #except:
        if 0:
            print(m)
            sys.exit(1)
        #print rmsd
        t += mrna.fn + ',' + str(rmsd) + '\n'
        #break    
    print(t.strip())
    if args.result:
        with open(args.result, 'w') as f:
            f.write(t)

        if args.sort:
            import pandas as pd
            df = pd.read_csv(args.result)
            df = df.sort_values('rmsd')
            df.to_csv(args.result, index=False)
            print(df.to_string(index=False))
        
        
