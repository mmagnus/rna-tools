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



def _read_fasta_sequence(fpath):
    """Return the first sequence found in a FASTA file as an uppercase string."""
    sequence_lines = []
    with open(fpath, 'r') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sequence_lines:
                    break
                continue
            sequence_lines.append(line)
    return ''.join(sequence_lines).upper()


def _map_fasta_onto_target_sequence(target_sequence, fasta_sequence, gap_open=None, gap_extend=None):
    """Align a FASTA sequence to the target sequence and return a residue-wise mapping."""
    if not target_sequence:
        raise ValueError('target structure does not contain a polymer sequence to align against')
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    gap_open_penalty = -abs(gap_open) if gap_open is not None else -1.0
    if gap_extend is None:
        gap_extend = gap_open
    gap_extend_penalty = -abs(gap_extend) if gap_extend is not None else -0.5
    aligner.open_gap_score = gap_open_penalty
    aligner.extend_gap_score = gap_extend_penalty

    alignments = aligner.align(target_sequence, fasta_sequence)
    try:
        alignment = alignments[0]
    except IndexError:
        raise ValueError('global alignment between target sequence and FASTA failed')

    seq1_line, _, seq2_line, _, _, _ = _build_alignment_strings(
        alignment, target_sequence, fasta_sequence)

    mapped = []
    for char_target, char_fasta in zip(seq1_line, seq2_line):
        if char_target == '-':
            # FASTA insertion w.r.t. structure; skip since no residue exists
            continue
        if char_fasta == '-':
            raise ValueError('FASTA sequence is missing residues present in the target structure')
        mapped.append(char_fasta)

    if len(mapped) != len(target_sequence):
        raise ValueError('alignment did not cover all target residues ({} mapped, {} expected)'.format(
            len(mapped), len(target_sequence)))

    return ''.join(mapped)


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

    try:
        seq1_ranges = alignment.aligned[0]
        seq2_ranges = alignment.aligned[1]
    except AttributeError:
        seq1_ranges = []
        seq2_ranges = []

    if len(seq1_ranges):
        seq1_start = int(seq1_ranges[0][0])
        seq1_end_idx = int(seq1_ranges[-1][1])
    else:
        seq1_start = seq1_idx
        seq1_end_idx = seq1_idx
    if len(seq2_ranges):
        seq2_start = int(seq2_ranges[0][0])
        seq2_end_idx = int(seq2_ranges[-1][1])
    else:
        seq2_start = seq2_idx
        seq2_end_idx = seq2_idx

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
    alignment_span = {
        'seq1_start': seq1_start,
        'seq2_start': seq2_start,
        'seq1_end': seq1_end_idx,
        'seq2_end': seq2_end_idx
    }
    return seq1_line, match_line, seq2_line, matched_pairs, matches, alignment_span


def _merge_alignment_columns(reference, existing_models, new_target, new_model, model_name='model'):
    """Merge pairwise alignments so every model shares the same gapped target string."""
    new_reference = []
    updated_models = [[] for _ in existing_models]
    aligned_new_model = []
    i = 0
    j = 0
    ref_len = len(reference)
    new_len = len(new_target)

    while i < ref_len or j < new_len:
        ref_char = reference[i] if i < ref_len else None
        new_char = new_target[j] if j < new_len else None

        if ref_char is None:
            new_reference.append(new_char)
            for existing in updated_models:
                existing.append('-')
            aligned_new_model.append(new_model[j])
            j += 1
            continue
        if new_char is None:
            new_reference.append(ref_char)
            for idx, existing in enumerate(updated_models):
                existing.append(existing_models[idx][i])
            aligned_new_model.append('-')
            i += 1
            continue

        if ref_char == new_char:
            new_reference.append(ref_char)
            for idx, existing in enumerate(updated_models):
                existing.append(existing_models[idx][i])
            aligned_new_model.append(new_model[j])
            i += 1
            j += 1
            continue

        if ref_char == '-':
            new_reference.append('-')
            for idx, existing in enumerate(updated_models):
                existing.append(existing_models[idx][i])
            aligned_new_model.append('-')
            i += 1
            continue

        if new_char == '-':
            new_reference.append('-')
            for existing in updated_models:
                existing.append('-')
            aligned_new_model.append(new_model[j])
            j += 1
            continue

        # mismatch: keep both columns but warn the user
        print('# warning: inconsistent target alignment between reference and {}: {!r} vs {!r}'.format(
            model_name, ref_char, new_char), file=sys.stderr)
        new_reference.append(ref_char)
        for idx, existing in enumerate(updated_models):
            existing.append(existing_models[idx][i])
        aligned_new_model.append('-')
        i += 1
        new_reference.append(new_char)
        for existing in updated_models:
            existing.append('-')
        aligned_new_model.append(new_model[j])
        j += 1
        continue

    merged_models = updated_models[:]
    merged_models.append(aligned_new_model)
    return new_reference, merged_models


def _print_alignment(seq1_line, match_line, seq2_line, header,
                     target_len=None, model_len=None, matches=None, residue_pairs=None,
                     rmsd=None, target_start_label=None, target_end_label=None,
                     model_start_label=None, model_end_label=None):
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
    if rmsd is not None:
        details.append('rmsd={}'.format(rmsd))
    if target_start_label:
        details.append('target_start={}'.format(target_start_label))
    if target_end_label:
        details.append('target_end={}'.format(target_end_label))
    if model_start_label:
        details.append('model_start={}'.format(model_start_label))
    if model_end_label:
        details.append('model_end={}'.format(model_end_label))
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

    def _sequence_position_label(self, index):
        """Return chain:resseq label for a given sequence index or None."""
        if index is None:
            return None
        try:
            idx = int(index)
        except (TypeError, ValueError):
            return None
        if idx < 0 or idx >= len(self.sequence_positions):
            return None
        return self.sequence_positions[idx]

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

        seq1_line, match_line, seq2_line, matched_pairs, matches, alignment_span = _build_alignment_strings(
            alignment, self.sequence, other_rnamodel.sequence)

        alignment_report = {
            'header': '# alignment between {} and {}'.format(self.fn, other_rnamodel.fn),
            'seq1_line': seq1_line,
            'match_line': match_line,
            'seq2_line': seq2_line,
            'target_len': len(seq1_line),
            'model_len': len(seq2_line),
            'matches': matches,
            'residue_pairs': len(matched_pairs),
            'seq1_start': alignment_span['seq1_start'],
            'seq1_end': alignment_span['seq1_end'],
            'seq2_start': alignment_span['seq2_start'],
            'seq2_end': alignment_span['seq2_end']
        }

        seq1_start_idx = alignment_span['seq1_start']
        seq1_end_idx = alignment_span['seq1_end']
        seq2_start_idx = alignment_span['seq2_start']
        seq2_end_idx = alignment_span['seq2_end']
        alignment_report['target_start_position'] = self._sequence_position_label(seq1_start_idx)
        alignment_report['target_end_position'] = self._sequence_position_label((seq1_end_idx - 1) if seq1_end_idx else None)
        alignment_report['model_start_position'] = other_rnamodel._sequence_position_label(seq2_start_idx)
        alignment_report['model_end_position'] = other_rnamodel._sequence_position_label((seq2_end_idx - 1) if seq2_end_idx else None)

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
            'alignment_pairs': len(matched_pairs),
            'alignment_report': alignment_report
        }
        return atoms_self, atoms_other, info

    def get_rmsd_to(self, other_rnamodel, way="", triple_mode=False, save=True):
        """Calc rmsd P-atom based rmsd to other rna model"""
        sup = Bio.PDB.Superimposer()
        alignment_report = None
        if args.align_sequence:
            if triple_mode:
                raise ValueError('Sequence alignment mode is not supported together with triple mode')
            self.atoms_for_rmsd, other_atoms_for_rmsd, info = self._atoms_from_sequence_alignment(other_rnamodel, way)
            alignment_report = info.get('alignment_report')
            other_rnamodel.alignment_report = alignment_report
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

            if alignment_report and getattr(args, 'print_alignment', False):
                _print_alignment(
                    alignment_report['seq1_line'],
                    alignment_report['match_line'],
                    alignment_report['seq2_line'],
                    alignment_report['header'],
                    target_len=alignment_report['target_len'],
                    model_len=alignment_report['model_len'],
                    matches=alignment_report['matches'],
                    residue_pairs=alignment_report['residue_pairs'],
                    rmsd=rms,
                    target_start_label=alignment_report.get('target_start_position'),
                    target_end_label=alignment_report.get('target_end_position'),
                    model_start_label=alignment_report.get('model_start_position'),
                    model_end_label=alignment_report.get('model_end_position'))

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


def expand_input_patterns(patterns):
    expanded = []
    for pattern in patterns:
        if any(char in pattern for char in '*?['):
            matches = glob.glob(pattern)
            if matches:
                expand_list = matches[:]
                sort_nicely(expand_list)
                expanded.extend(expand_list)
            else:
                expanded.append(pattern)
        else:
            expanded.append(pattern)
    return expanded


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t',"--target",
                         help="target file")
    parser.add_argument('--target-fasta',
                        help='FASTA file with the sequence that should override the target PDB sequence for alignments')
    parser.add_argument('--early-stop', type=int, metavar='N', default=None,
                        help='stop after processing N models (useful for debugging); default: process all models')
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
    parser.add_argument('--alignment-fasta', nargs='?', const='seq.fasta',
                        help='write the aligned target/model sequences to FASTA; optional output path, default: seq.fasta')
    parser.add_argument('--add-rmsd-to-fasta-header', action='store_true',
                        help='prefix each FASTA model header with its RMSD value when writing --alignment-fasta')
    parser.add_argument('--sort-by-rmsd', action='store_true',
                        help='order FASTA entries by descending RMSD when using --alignment-fasta')
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
    if args.alignment_fasta and not args.align_sequence:
        parser.error('--alignment-fasta requires --align-sequence')

    target = RNAmodel(args.target)
    if args.target_fasta:
        fasta_sequence = _read_fasta_sequence(args.target_fasta)
        if not fasta_sequence:
            parser.error('target FASTA does not contain any sequence data')
        polymer_len = len(getattr(target, 'sequence_atom_maps', []))
        if polymer_len == 0:
            parser.error('target FASTA override requires the target structure to contain polymer residues')
        try:
            mapped_sequence = _map_fasta_onto_target_sequence(
                target.sequence,
                fasta_sequence,
                args.align_gap_open,
                args.align_gap_extend
            )
        except ValueError as exc:
            parser.error('target FASTA override failed: {}'.format(exc))
        if len(mapped_sequence) != polymer_len:
            parser.error('internal error: remapped FASTA length ({}) does not match target residues ({})'.format(
                len(mapped_sequence), polymer_len))
        print('# overriding target sequence with FASTA {} (length={} -> mapped {})'.format(
            args.target_fasta, len(fasta_sequence), len(mapped_sequence)), file=sys.stderr)
        target.sequence = mapped_sequence

    target_sequence = target.sequence if target.sequence else ''
    target_len = len(target_sequence)

    if args.print_target_sequence:
        print('# target sequence (polymer residues only):', target.fn)
        if target.sequence:
            print(target.sequence)
            print('# length:', len(target.sequence))
        else:
            print('# warning: target structure does not contain polymer residues with standard names')

    models = expand_input_patterns(args.files)
    tmp = []
    if args.ignore_files:
        for f in models:
            if args.debug: print(f)
            if args.ignore_files in f:
                continue
            tmp.append(f)
        models = tmp
 
    print('# of models:', len(models))
    processed = 0
    header_columns = ['fn', args.column_name]
    if args.align_sequence:
        header_columns.extend(['target_alignment_start', 'target_alignment_end'])
    t = ','.join(header_columns) + '\n'
    alignment_entries = [] if args.alignment_fasta else None
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
        alignment = getattr(mrna, 'alignment_report', None)
        #print rmsd
        if alignment_entries is not None:
            if alignment and alignment.get('seq1_line') and alignment.get('seq2_line'):
                alignment_entries.append({
                    'model_name': mrna.fn,
                    'target_seq': alignment['seq1_line'],
                    'model_seq': alignment['seq2_line'],
                    'seq1_start': alignment.get('seq1_start', 0),
                    'seq1_end': alignment.get('seq1_end', target_len),
                    'target_start_position': alignment.get('target_start_position'),
                    'target_end_position': alignment.get('target_end_position'),
                    'rmsd': rmsd
                })
            else:
                print('# warning: no alignment info collected for', mrna.fn, file=sys.stderr)
        row = [mrna.fn, str(rmsd)]
        if args.align_sequence:
            start_label = alignment.get('target_start_position') if alignment else ''
            end_label = alignment.get('target_end_position') if alignment else ''
            row.extend([
                start_label if start_label is not None else '',
                end_label if end_label is not None else ''
            ])
        t += ','.join(row) + '\n'
        processed += 1
        if args.early_stop and processed >= args.early_stop:
            print('# early stop after {} models'.format(processed), file=sys.stderr)
            break
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

    if alignment_entries is not None:
        fasta_path = args.alignment_fasta
        valid_entries = [entry for entry in alignment_entries if entry['target_seq'] and entry['model_seq']]
        if not valid_entries:
            print('# warning: requested FASTA output but no alignment data was generated', file=sys.stderr)
        else:
            padded_entries = []
            for entry in valid_entries:
                seq1_start = entry.get('seq1_start', 0)
                seq1_end = entry.get('seq1_end', target_len)
                target_leading = target_sequence[:seq1_start] if target_sequence else '-' * max(seq1_start, 0)
                target_trailing = target_sequence[seq1_end:] if target_sequence else '-' * max((target_len - seq1_end) if target_len else 0, 0)

                leading = '-' * len(target_leading)
                trailing = '-' * len(target_trailing)
                padded_entries.append({
                    'model_name': entry['model_name'],
                    'target_seq': target_leading + entry['target_seq'] + target_trailing,
                    'model_seq': leading + entry['model_seq'] + trailing,
                    'rmsd': entry.get('rmsd')
                })

            if args.sort_by_rmsd:
                def rmsd_key(entry):
                    value = entry.get('rmsd')
                    try:
                        return float(value)
                    except (TypeError, ValueError):
                        return float('-inf')
                padded_entries = sorted(padded_entries, key=rmsd_key, reverse=True)

            ref_target = list(padded_entries[0]['target_seq'])
            model_names = [padded_entries[0]['model_name']]
            model_columns = [list(padded_entries[0]['model_seq'])]

            for entry in padded_entries[1:]:
                ref_target, model_columns = _merge_alignment_columns(
                    ref_target,
                    model_columns,
                    list(entry['target_seq']),
                    list(entry['model_seq']),
                    entry['model_name']
                )
                model_names.append(entry['model_name'])

            with open(fasta_path, 'w') as fasta:
                fasta.write('>{}\n'.format(target.fn))
                fasta.write(''.join(ref_target) + '\n')
                for name, seq in zip(model_names, model_columns):
                    header_name = name
                    if args.add_rmsd_to_fasta_header:
                        rmsd_value = next((entry['rmsd'] for entry in padded_entries if entry['model_name'] == name), None)
                        if rmsd_value is not None:
                            header_name = '{}|rmsd={}'.format(rmsd_value, name)
                    fasta.write('>{}\n'.format(header_name))
                    fasta.write(''.join(seq) + '\n')
            print('# FASTA alignment written to', fasta_path)

        
def _read_fasta_sequence(fpath):
    """Return the first sequence found in a FASTA file as an uppercase string."""
    sequence_lines = []
    with open(fpath, 'r') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sequence_lines:
                    break
                continue
            sequence_lines.append(line)
    return ''.join(sequence_lines).upper()
