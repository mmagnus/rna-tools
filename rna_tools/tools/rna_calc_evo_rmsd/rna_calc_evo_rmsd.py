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

All vs all: with ``--all_vs_all`` every structure is mapped on its best Rfam
family, structures of each family are aligned together to the family covariance
model and RMSD and sequence identity (``esl-alipid`` of Easel) are calculated for
all pairs within each family. A plot of sequence identity vs RMSD is saved to
``<output>_seqid_vs_rmsd.png``::

    rna_calc_evo_rmsd.py --all_vs_all --rfam_db Rfam.cm -o rmsd.csv *.pdb

Sequence identity is calculated in all modes (column ``seq_identity``) with
``esl-alipid`` (Easel, https://github.com/EddyRivasLab/easel), if available.

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
pd.set_option('display.max_columns', None)
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
    parser.add_argument('-t', "--target", help="the native structure file (required, unless --all_vs_all)")
    parser.add_argument('--all_vs_all', action='store_true',
                        help="map all structures (files and --target, if given) on Rfam families, build an alignment per family \
                        (<output>_rfam/<family>.sto) and calculate RMSD and sequence identity for all pairs within each family; \
                        see the plot <output>_seqid_vs_rmsd.png")
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
_esl_alipid_warned = False


def calc_seq_identity(aligned1, aligned2, verbose=False):
    """Calculate sequence identity (%id) of two aligned sequences with esl-alipid (Easel).

    %id = identical aligned residues / min(length of seq1, length of seq2) * 100,
    as defined by esl-alipid. Returns NaN if esl-alipid is not available.

    esl-alipid is one of Easel miniapps (https://github.com/EddyRivasLab/easel),
    it is installed with HMMER/Infernal source distributions."""
    global _esl_alipid_warned
    if not shutil.which('esl-alipid'):
        if not _esl_alipid_warned:
            print(' WARNING: esl-alipid (Easel) not found, sequence identity not calculated')
            _esl_alipid_warned = True
        return float('nan')
    # write a clean two-seq Stockholm file, Easel is strict about the format
    def clean(seq):
        return ''.join('-' if c in '.~' else c for c in seq)
    with tempfile.NamedTemporaryFile('w', suffix='.sto', delete=False) as f:
        f.write('# STOCKHOLM 1.0\n\nseq1 %s\nseq2 %s\n//\n' % (clean(aligned1), clean(aligned2)))
        fn = f.name
    try:
        out = RfamAlign.run(['esl-alipid', '--rna', '--noheader', fn], verbose)
    finally:
        os.remove(fn)
    # seqname1 seqname2 %id nid denomid %match nmatch denommatch
    for line in out.splitlines():
        cols = line.split()
        if len(cols) >= 3 and cols[0] == 'seq1' and cols[1] == 'seq2':
            return float(cols[2])
    raise Exception('Could not parse esl-alipid output:\n' + out)


def align_structures(target_struc, model_struc, verbose=False):
    """Align sequences of two structures with a global pairwise alignment.

    End gaps are penalized less than internal gaps, so a fragment can be aligned to a
    full-length structure. This is a sequence-only alignment, for divergent homologs
    provide a (structural) alignment with -a instead.

    :returns: (target_positions, model_positions, (aligned target seq, aligned model seq)),
              1-based positions of the aligned (non-gap) residues, two lists of the same length"""
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
    return pos1, pos2, (aln[0], aln[1])


def rfam_align(targetfn, files, output_dir, rfam_db=None, rfam_evalue=None, verbose=False):
    """Align the target and models using Rfam, see RfamAlign.py.

    :returns: list of (target fn, model fn, target positions, model positions, family, (aligned target seq, aligned model seq))"""
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
                pairs.append((targetfn, fn, target_pos, model_pos, family, (ra.get_seq(target_name), ra.get_seq(name))))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return pairs


def rfam_all_vs_all(files, output_dir, rfam_db=None, rfam_evalue=None, verbose=False):
    """Map structures on Rfam families and align structures of each family, see RfamAlign.py.

    Each structure is assigned to its best Rfam family (the highest score), structures
    of each family are aligned together to the family CM (one alignment per family,
    saved to output_dir/<family>.sto), and all pairs within a family are returned.

    :returns: list of (target fn, model fn, target positions, model positions, family,
              (aligned target seq, aligned model seq))"""
    rfam_db = RfamAlign.get_rfam_db(rfam_db)
    print(' Rfam db:', rfam_db)
    names = RfamAlign.get_seq_names(files)
    seqs = [(name, get_sequence(RNAmodel.parse(fn))) for name, fn in zip(names, files)]
    seq_of = dict(seqs)

    pairs = []
    workdir = tempfile.mkdtemp()
    try:
        hits = RfamAlign.cmscan(seqs, rfam_db, workdir, rfam_evalue, verbose)
        scores = RfamAlign.get_best_scores(hits)
        families = {}  # family: [(fn, name), ...]
        for fn, name in zip(files, names):
            fams = scores.get(name)
            if not fams:
                print(' WARNING: no Rfam hit, skipped:', fn)
                continue
            family = max(fams, key=fams.get)
            print(' %s: %s' % (name, family))
            families.setdefault(family, []).append((fn, name))

        if families and not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        for family, structures in sorted(families.items()):
            if len(structures) < 2:
                print(' WARNING: only one structure in %s, skipped: %s' % (family, structures[0][0]))
                continue
            aln_fn = os.path.join(output_dir, family + '.sto')
            RfamAlign.cmalign(family, [(name, seq_of[name]) for fn, name in structures],
                              rfam_db, workdir, aln_fn, verbose)
            print(' alignment saved:', aln_fn)
            ra = RNAalignment(aln_fn, verbose=verbose)
            for i, (fn1, name1) in enumerate(structures):
                for fn2, name2 in structures[i + 1:]:
                    pos1, pos2 = ra.get_paired_positions(name1, name2, verbose=verbose)
                    pairs.append((fn1, fn2, pos1, pos2, family, (ra.get_seq(name1), ra.get_seq(name2))))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return pairs


def get_rmsd_table(pairs, group_name='', verbose=False):
    """Calculate RMSD and sequence identity for pairs.

    :param pairs: list of (target fn, model fn, target positions, model positions, family,
                  (aligned target seq, aligned model seq))
    :returns: pandas DataFrame"""
    data = {'target': [], 'model': [], 'rmsd': [], 'n_residues': [], 'seq_identity': [], 'family': [], 'group_name': []}
    for targetfn, f, target_pos, model_pos, family, aligned in pairs:
        target = RNAmodel(targetfn, target_pos, save=False, output_dir=None)
        model = RNAmodel(f, model_pos, save=False, output_dir=None)
        rmsd = target.get_rmsd_to(model)
        data['target'].append(target)
        data['model'].append(model)
        data['rmsd'].append(rmsd)
        data['n_residues'].append(len(target_pos))
        data['seq_identity'].append(calc_seq_identity(aligned[0], aligned[1], verbose))
        data['family'].append(family)
        data['group_name'].append(group_name)
    return pd.DataFrame(data, columns=('target', 'model', 'rmsd', 'n_residues', 'seq_identity', 'family', 'group_name'))


def plot_seq_identity_vs_rmsd(df, fn):
    """Scatter plot of sequence identity vs RMSD, one color per Rfam family, save to fn."""
    df = df.dropna(subset=['seq_identity'])
    if df.empty:
        return
    fig, ax = plt.subplots()
    for family, g in df.groupby('family'):
        ax.scatter(g['seq_identity'], g['rmsd'], label=family or None, alpha=0.8)
    if (df['family'] != '').any():
        ax.legend(title='Rfam family')
    ax.set_xlabel('sequence identity [%] (esl-alipid)')
    ax.set_ylabel('RMSD [A]')
    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0)
    fig.savefig(fn, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(' plot saved:', fn)


def save_results(df, output_fn):
    """Save the table (csv), the plot of sorted RMSDs (.png) and of sequence identity vs RMSD (_seqid_vs_rmsd.png)."""
    df.to_csv(output_fn)
    base = os.path.splitext(output_fn)[0]
    try:
        df.sort_values(by='rmsd').plot(y='rmsd', use_index=False)
    except TypeError:
        print('Check if the representives have tags, e.g. c1_thf_pk_...')
    plt.savefig(base + '.png')
    plt.close('all')
    plot_seq_identity_vs_rmsd(df, base + '_seqid_vs_rmsd.png')


def calc_evo_rmsd_all_vs_all(files, output_fn=None, group_name='', verbose=False, rfam_db=None, rfam_evalue=None):
    """Map structures on Rfam families, align them and calculate RMSD and sequence
    identity for all pairs of structures within each family (see rfam_all_vs_all)."""
    global _atom_summary_printed
    if not _atom_summary_printed:
        print(get_atom_selection_summary())
        _atom_summary_printed = True
    files = list(dict.fromkeys(files))  # remove duplicates, keep order
    output_dir = os.path.splitext(output_fn or 'evoclust_rmsd.csv')[0] + '_rfam'
    pairs = rfam_all_vs_all(files, output_dir, rfam_db, rfam_evalue, verbose)
    df = get_rmsd_table(pairs, group_name, verbose)
    if output_fn:
        save_results(df, output_fn)
    return df


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

    pairs = []  # see get_rmsd_table
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
                target_pos, model_pos, aligned = align_structures(target_struc, RNAmodel.parse(f), verbose)
                if not target_pos:
                    raise Exception('Sequences of %s and %s could not be aligned' % (targetfn, f))
                pairs.append((targetfn, f, target_pos, model_pos, '', aligned))
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
                    pairs.append((targetfn, f, target_pos, model_pos, '',
                                  (ra.get_seq(target_name_alignment), ra.get_seq(rs_name_alignment))))

    df = get_rmsd_table(pairs, group_name, verbose)
    if output_fn:
        save_results(df, output_fn)
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
    if opts.all_vs_all:
        try:
            df = calc_evo_rmsd_all_vs_all(([opts.target] if opts.target else []) + opts.files, opts.output_fn,
                                          opts.group_name, opts.verbose, opts.rfam_db, opts.rfam_evalue)
        except RfamAlign.RfamAlignError as e:
            print('Error:', e)
            sys.exit(1)
        print(df)
        sys.exit(0)
    if not opts.target:
        parser.error('the following arguments are required: -t/--target (or use --all_vs_all)')
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
