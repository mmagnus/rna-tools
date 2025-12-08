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
from RNAmodel import RNAmodel, get_atom_selection_summary
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
                        help="Stockholm alignment file with either an x/EvoClust selector line, a #=GC RF reference annotation, or gap-free columns that will be auto-detected (e.g. test_data/rp14sub.stk)",
                        required=True)
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


def calc_evo_rmsd(targetfn, target_name_alignment, files, mapping_fn, rna_alignment_fn, group_name='', output_fn=None,
                  verbose=False):
    global _atom_summary_printed
    if not _atom_summary_printed:
        print(get_atom_selection_summary())
        _atom_summary_printed = True
    ra = RNAalignment(rna_alignment_fn, verbose=verbose)
    print('target', targetfn)
    target_residues = ra.get_range(target_name_alignment, verbose=verbose)
    target = RNAmodel(targetfn, target_residues, save=False, output_dir=None)

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
    models = []

    for rs in rnastruc:
        try:
            rs_name_alignment, rs_name_dir = [value.strip() for value in rs.split(':', 1)]  # target:rp14_farna_eloop_nol2fixed_cst
        except ValueError:
            # if -m 'tpp|tpp_pdb|CP000050.1/ ..
            # rnastruc: ['tpp', 'tpp_pdb', 'CP000050.1/1019813-1019911:tc5_pdb', 'AE017180.1/640928-641029:tae_pdb', 'BX248356.1/234808-234920:tb2_pdb']
            # Traceback (most recent call last):
            # File "/home/magnus/work/src/evoClustRNA/evoClustRNA.py", line 100, in <module>
            # raise Exception("There is an error in your mapping, check all : and | carefully")
            # Exception: There is an error in your mapping, check all : and | carefully
            raise Exception("There is an error in your mapping, check all : and | carefully")

        # print ' ', rs_name_alignment,'<->', rs_name_dir # AACY023581040 <-> aacy23_cst
        for f in files:
            #
            if rs_name_dir in f:  # rp14_farna_eloop_nol2fixed_cst*pdb
                # print(rs_name_dir)
                # models.extend(get_rna_models_from_dir(input_dir + os.sep + rs_name_dir, ra.get_range(rs_name_alignment), False, False)[:])
                # print(ra.get_range(rs_name_alignment))
                models.append(RNAmodel(f, ra.get_range(rs_name_alignment), False, False))

    data = {'target': [], 'model': [], 'rmsd': [], 'group_name': []}
    for model in models:
        # print
        # print r1.fn, r2.fn, r1.get_rmsd_to(r2)#, 'tmp.pdb')
        rmsd = target.get_rmsd_to(model)  # , 'tmp.pdb')
        # print target, model, rmsd, group_name
        data['target'].append(target)
        data['model'].append(model)
        data['rmsd'].append(rmsd)
        data['group_name'].append(group_name)
    df = pd.DataFrame(data, columns=('target', 'model', 'rmsd', 'group_name'))
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
    if not target_name:
        target_name = os.path.splitext(os.path.basename(opts.target))[0]
        print(' target name not provided; using basename:', target_name)
    df = calc_evo_rmsd(opts.target, target_name, opts.files, opts.mapping_fn,
                       opts.rna_alignment_fn, opts.group_name, opts.output_fn,
                       verbose=opts.verbose)
    print(df)
