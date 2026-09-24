#!/usr/bin/env python

"""Calculate RMSD between structures based on a Stockholm alignment and selector lines.

When RNA models are loaded, models ending with 'template.pdb' are ignore.

c1_tha_96cdea07- tha in mapping

Structure discovery
-------------------

The script accepts individual PDB paths as positional arguments.  When a
mapping file is supplied through ``--mapping_fn`` every entry is expected to be
``<alignment_id>:<substring>`` and the substring is matched against the provided
PDB file paths to decide which structures belong to a given alignment sequence.
If the mapping file is omitted the script automatically builds
``<basename>:<basename>`` pairs from the positional PDB filenames, effectively
using the PDB basename both as the alignment identifier and as the lookup
substring.

The alignment needs either an explicit ``x``/``EvoClust`` sequence, a
``#=GC RF`` reference annotation, or columns that are gap-free across all
sequences (an x-line will be inferred) to indicate which positions should be
used for the RMSD measurement.

When ``--target_name`` is not given the script uses the basename of the target
structure passed with ``--target`` as the identifier in the alignment.

Automated mode (no alignment), Rfam
-----------------------------------

If ``-a/--rna_alignment_fn`` is omitted, sequences are taken from the PDB files
and searched against Rfam with ``cmscan``. For each model, the Rfam family hit
by both the target and the model (the best one, if there are more) is taken,
the sequences are aligned to its covariance model with ``cmalign``, and the
residues in the consensus (``#=GC RF``) columns are used for the RMSD
calculation. No alignment, mapping or selector line has to be prepared by hand::

    rna_calc_evo_rmsd.py --rfam_db Rfam.cm -t test_data/1ehz_std.pdb test_data/6Y2L_2_std.pdb

The alignments are saved to ``<output>_rfam/<family>.sto`` (they can be
inspected and re-used with ``-a``). Models that do not share a family with the
target are skipped. Requires Infernal and Rfam.cm (see RfamAlign.py).

With ``--auto pairwise`` the sequences are instead aligned with a simple
pairwise sequence alignment (no Rfam needed), and all aligned positions are
used. This works best for closely related sequences.

Residues are always paired according to the alignment columns (the n-th residue
in the alignment is the n-th nucleotide in the PDB file), not according to the
residue numbers in the PDB files.
"""
from __future__ import print_function
import pandas as pd
pd.set_option('display.width', 1000)
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import argparse

import sys
import math
import glob
import re
import os
import shutil

import Bio.PDB.PDBParser
import Bio.PDB.Superimposer
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO, Superimposer

from RNAalignment import RNAalignment
from RNAmodel import RNAmodel, get_atom_selection_summary, get_sequence
from Bio.Align import PairwiseAligner
import RfamAlign
import tempfile
import csv

debug = False


def get_rna_models_from_dir(directory, residues, save, output_dir):
    models = []
    if not os.path.exists(directory):
        raise Exception('Dir does not exist! ', directory)
    files = glob.glob(directory + "/*.pdb")
    files_sorted = sort_nicely(files)
    for f in files_sorted:
        # ignore files that can be found in your folder
        # be careful with this --magnus
        if f.endswith('template.pdb'):
            continue
        if 'clust01X' in f:
            continue
        if 'clust02X' in f:
            continue
        if 'clust03X' in f:
            continue
        models.append(RNAmodel(f, residues, save, output_dir))
    return models


def sort_nicely(l):
    """ Sort the given list in the way that humans expect.

    http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
    """
    def convert(text): return int(text) if text.isdigit() else text

    def alphanum_key(key): return [convert(c) for c in re.split('([0-9]+)', key)]
    l.sort(key=alphanum_key)
    return l


def parse_num_list(s):
    """ http://stackoverflow.com/questions/6512280/accept-a-range-of-numbers-in-the-form-of-0-5-using-pythons-argparse """
    m = re.match(r'(\d+)(?:-(\d+))?$', str(s))
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        return s
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start, 10), int(end, 10) + 1))

# def pair:
#    def __init__(self, f1, f2):
#        self.f1
#        self.f2
#    def calc_distance():
#        pass


def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', "--rna_alignment_fn",
                        help="Stockholm alignment file with either an x/EvoClust selector line, a #=GC RF reference annotation, or gap-free columns that will be auto-detected (e.g. test_data/rp14sub.stk). \
                        If omitted, the alignment is made automatically, see --auto.",
                        default=None)
    parser.add_argument('--auto', choices=['rfam', 'pairwise'], default='rfam',
                        help="how to align the structures if no alignment (-a) is given: rfam - find the Rfam family shared by the target and a model (cmscan) and align them to its covariance model (cmalign), \
                        consensus (RF) columns are used; pairwise - a pairwise sequence alignment of the model to the target, all aligned residues are used (default: rfam)")
    parser.add_argument('--rfam_db', help="path to Rfam.cm (cmpress'ed), by default $RFAM_DB_PATH or RFAM_DB_PATH from the rna-tools config")
    parser.add_argument('--rfam_evalue', type=float,
                        help="report Rfam hits with E-value <= this value; by default Rfam gathering thresholds (--cut_ga) are used")
    parser.add_argument('-t', "--target", help="the native structure file", required=True)
    parser.add_argument('-o', "--output_fn", help="output csv file", default="evoclust_rmsd.csv")
    parser.add_argument('-n', "--target_name",
                        help="target name in the alignment, used to map target on the alignment, e.g. target, ade, rp14 etc. Defaults to the basename of --target.")
    parser.add_argument('-m', "--mapping_fn", help="map folders on the drive with sequence names in the alignment (<name in the alignment>:<folder name>), use | to \
    for multiple seqs, e.g. 'target:rp14_farna_eloop_nol2fixed_cst|AACY023581040:aacy23_cst', use | as a separator. If omitted, PDB basenames are used.",
                        default=None)
    parser.add_argument('files', nargs='+', help='files')
    parser.add_argument('-g', '--group_name',
                        help='name given group of structure, helps to analyze results', default='')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print alignment records and the selector (x-line) used for RMSD calculations')
    return parser

_atom_summary_printed = False


def align_structures(target_struc, model_struc, verbose=False):
    """Align sequences of two structures with a global pairwise alignment.

    End gaps are penalized less than internal gaps, so a fragment can be aligned to a
    full-length structure. This is a sequence-only alignment, for divergent homologs
    provide a (structural) alignment with -a instead.

    :returns: (target_positions, model_positions), 1-based positions of the aligned
              (non-gap) residues, two lists of the same length"""
    seq1 = get_sequence(target_struc)
    seq2 = get_sequence(model_struc)
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aligner.end_gap_score = -1
    aln = aligner.align(seq1, seq2)[0]
    pos1, pos2 = [], []
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        pos1.extend(range(s1 + 1, e1 + 1))
        pos2.extend(range(s2 + 1, e2 + 1))
    if verbose:
        print(aln)
    return pos1, pos2


def rfam_align(targetfn, files, output_dir, rfam_db=None, rfam_evalue=None, verbose=False):
    """Align the target and models using Rfam, see RfamAlign.py.

    :returns: list of (model fn, target positions, model positions, family)"""
    rfam_db = RfamAlign.get_rfam_db(rfam_db)
    print(' Rfam db:', rfam_db)
    fns = [targetfn] + files
    names = RfamAlign.get_seq_names(fns)
    target_name = names[0]
    seqs = [(name, get_sequence(RNAmodel.parse(fn))) for name, fn in zip(names, fns)]

    pairs = []
    workdir = tempfile.mkdtemp()
    try:
        hits = RfamAlign.cmscan(seqs, rfam_db, workdir, rfam_evalue, verbose)
        scores = RfamAlign.get_best_scores(hits)
        print(' Rfam families of the target (%s): %s' % (target_name, ', '.join(sorted(scores.get(target_name, {}))) or 'none'))
        if not scores.get(target_name):
            raise RfamAlign.RfamAlignError('No Rfam hit for the target %s' % targetfn)
        families = {}  # family: [(fn, name), ...]
        for fn, name in zip(fns[1:], names[1:]):
            family = RfamAlign.find_common_family(scores, target_name, name)
            if not family:
                print(' WARNING: no Rfam family shared with the target, skipped:', fn)
                continue
            print(' %s: %s' % (name, family))
            families.setdefault(family, []).append((fn, name))

        if families and not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        seq_of = dict(seqs)
        for family, models in families.items():
            aln_fn = os.path.join(output_dir, family + '.sto')
            RfamAlign.cmalign(family, [(target_name, seq_of[target_name])] + [(name, seq_of[name]) for fn, name in models],
                              rfam_db, workdir, aln_fn, verbose)
            print(' alignment saved:', aln_fn)
            ra = RNAalignment(aln_fn, verbose=verbose)
            for fn, name in models:
                target_pos, model_pos = ra.get_paired_positions(target_name, name, verbose=verbose)
                pairs.append((fn, target_pos, model_pos, family))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return pairs


def calc_evo_rmsd(targetfn, target_name_alignment, files, mapping_fn, rna_alignment_fn, group_name='', output_fn=None,
                  verbose=False, auto='rfam', rfam_db=None, rfam_evalue=None):
    """Calculate RMSD of models to the target.

    If rna_alignment_fn is None, the models are automatically aligned to the target
    and mapping_fn is ignored: auto='rfam' (see rfam_align) or auto='pairwise'
    (see align_structures)."""
    global _atom_summary_printed
    if not _atom_summary_printed:
        print(get_atom_selection_summary())
        _atom_summary_printed = True
    print('target', targetfn)

    pairs = []  # (model fn, target positions, model positions, family)
    if not rna_alignment_fn:
        files = [f for f in files if os.path.abspath(f) != os.path.abspath(targetfn)]
        if auto == 'rfam':
            print(' alignment not provided; aligning models to the target with Rfam')
            output_dir = os.path.splitext(output_fn or 'evoclust_rmsd.csv')[0] + '_rfam'
            pairs = rfam_align(targetfn, files, output_dir, rfam_db, rfam_evalue, verbose)
        else:
            print(' alignment not provided; aligning models to the target (pairwise)')
            target_struc = RNAmodel.parse(targetfn)
            for f in files:
                if verbose:
                    print(' ', os.path.basename(targetfn), '<->', os.path.basename(f))
                target_pos, model_pos = align_structures(target_struc, RNAmodel.parse(f), verbose)
                if not target_pos:
                    raise Exception('Sequences of %s and %s could not be aligned' % (targetfn, f))
                pairs.append((f, target_pos, model_pos, ''))
    else:
        ra = RNAalignment(rna_alignment_fn, verbose=verbose)
        # check if the target is in the alignment
        ra.get_range(target_name_alignment, verbose=verbose)

        # parse mapping to get models (list models)
        if mapping_fn:
            mapping_content = open(mapping_fn).read().replace('\n', '').strip()
            rnastruc = [item.strip() for item in mapping_content.split('|') if item.strip()]
        else:
            rnastruc = []
            for pdb_path in files:
                pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
                rnastruc.append(f"{pdb_name}:{pdb_name}")
            print(' mapping file not provided; using PDB basenames as alignment IDs')
        print(' # of rnastruc :', len(rnastruc))
        print(' rnastruc:', rnastruc)
        print(' WARNING: if any of your PDB file is missing, check mapping!')

        for rs in rnastruc:
            try:
                rs_name_alignment, rs_name_dir = [value.strip() for value in rs.split(':', 1)]  # target:rp14_farna_eloop_nol2fixed_cst
            except ValueError:
                # if -m 'tpp|tpp_pdb|CP000050.1/ ..
                # rnastruc: ['tpp', 'tpp_pdb', 'CP000050.1/1019813-1019911:tc5_pdb', 'AE017180.1/640928-641029:tae_pdb', 'BX248356.1/234808-234920:tb2_pdb']
                raise Exception("There is an error in your mapping, check all : and | carefully")

            # print ' ', rs_name_alignment,'<->', rs_name_dir # AACY023581040 <-> aacy23_cst
            for f in files:
                if rs_name_dir in f:  # rp14_farna_eloop_nol2fixed_cst*pdb
                    target_pos, model_pos = ra.get_paired_positions(target_name_alignment, rs_name_alignment,
                                                                    verbose=verbose)
                    pairs.append((f, target_pos, model_pos, ''))

    data = {'target': [], 'model': [], 'rmsd': [], 'n_residues': [], 'family': [], 'group_name': []}
    for f, target_pos, model_pos, family in pairs:
        target = RNAmodel(targetfn, target_pos, save=False, output_dir=None)
        model = RNAmodel(f, model_pos, save=False, output_dir=None)
        rmsd = target.get_rmsd_to(model)
        data['target'].append(target)
        data['model'].append(model)
        data['rmsd'].append(rmsd)
        data['n_residues'].append(len(target_pos))
        data['family'].append(family)
        data['group_name'].append(group_name)
    df = pd.DataFrame(data, columns=('target', 'model', 'rmsd', 'n_residues', 'family', 'group_name'))
    if output_fn:
        df.to_csv(output_fn)
        try:
            df.sort_values(by='rmsd').plot(y='rmsd', use_index=False)
        except TypeError:
            print('Check if the representives have tags, e.g. c1_thf_pk_...')
        plt.savefig(output_fn.replace('.csv', '.png'))
    return df


def test():
    mapping = 'target:rp14_farna_eloop_nol2fixed_cst|AACY023581040:aacy23_cst|AJ630128:aj63_cst'
    x = calc_evo_rmsd("test_data/rp14/rp14_5ddp_bound_clean_ligand.pdb", 'target',
                      ['test_data/rp14/rp14_farna_eloop_nol2fixed_cst/rp14_farna_eloop_nol2fixed_cst.out.1.pdb'],
                      mapping, rna_alignment_fn="test_data/rp14/rp14sub.stk")
    print(x)
    sys.exit(0)


# main
if __name__ == '__main__':
    # if True: test()
    parser = get_parser()
    opts = parser.parse_args()
    target_name = opts.target_name
    if not target_name and opts.rna_alignment_fn:
        target_name = os.path.splitext(os.path.basename(opts.target))[0]
        print(' target name not provided; using basename:', target_name)
    try:
        df = calc_evo_rmsd(opts.target, target_name, opts.files, opts.mapping_fn,
                           opts.rna_alignment_fn, opts.group_name, opts.output_fn,
                           verbose=opts.verbose, auto=opts.auto, rfam_db=opts.rfam_db,
                           rfam_evalue=opts.rfam_evalue)
    except RfamAlign.RfamAlignError as e:
        print('Error:', e)
        sys.exit(1)
    print(df)
