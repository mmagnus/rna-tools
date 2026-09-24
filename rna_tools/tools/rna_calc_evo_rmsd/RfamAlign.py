#!/usr/bin/env python
"""Find the Rfam family shared by structures and align them to its covariance model.

1. ``cmscan`` all sequences against Rfam (Rfam.cm),
2. for each model, pick the best family hit by both the target and the model,
3. ``cmfetch`` the family CM and ``cmalign`` the sequences to it.

The resulting Stockholm alignment has a ``#=GC RF`` line, so the consensus
(RF) columns are used to pair residues for the RMSD calculation.

Setup:

- install Infernal http://eddylab.org/infernal/ (cmscan, cmfetch, cmalign in $PATH),
- download ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz, gunzip and run ``cmpress Rfam.cm``,
- give the path with ``--rfam_db``, or set ``RFAM_DB_PATH`` (environment variable or the rna-tools config).
"""
from __future__ import print_function
import os
import re
import shutil
import subprocess


class RfamAlignError(Exception):
    pass


def get_rfam_db(rfam_db=None):
    """Get path to Rfam.cm: given path, $RFAM_DB_PATH or RFAM_DB_PATH from rna_tools_config."""
    if not rfam_db:
        rfam_db = os.environ.get('RFAM_DB_PATH')
    if not rfam_db:
        try:
            from rna_tools.rna_tools_config import RFAM_DB_PATH
            rfam_db = RFAM_DB_PATH
        except ImportError:
            pass
    if not rfam_db or not os.path.isfile(rfam_db):
        raise RfamAlignError('Rfam.cm not found (%s). Give it with --rfam_db or set RFAM_DB_PATH. '
                             'Download: ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz, '
                             'gunzip and run cmpress Rfam.cm' % rfam_db)
    for tool in ('cmscan', 'cmfetch', 'cmalign'):
        if not shutil.which(tool):
            raise RfamAlignError('%s not found, install Infernal http://eddylab.org/infernal/' % tool)
    return rfam_db


def run(cmd, verbose=False):
    if verbose:
        print(' ' + ' '.join(cmd))
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if p.returncode != 0:
        raise RfamAlignError('%s failed:\n%s%s' % (cmd[0], p.stdout, p.stderr))
    return p.stdout


def write_fasta(seqs, fn):
    with open(fn, 'w') as f:
        for name, seq in seqs:
            f.write('>%s\n%s\n' % (name, seq))


def cmscan(seqs, rfam_db, workdir, evalue=None, verbose=False):
    """Run cmscan of seqs against Rfam.

    :param seqs: list of (name, seq)
    :param evalue: if None, use Rfam gathering thresholds (--cut_ga), otherwise report hits with E-value <= evalue
    :returns: list of hits, dicts with keys: family, accession, name, seq_from, seq_to, strand, score, evalue
    """
    fasta = os.path.join(workdir, 'seqs.fa')
    tblout = os.path.join(workdir, 'cmscan.tblout')
    write_fasta(seqs, fasta)
    cmd = ['cmscan', '--rfam', '--nohmmonly', '--tblout', tblout]
    if evalue is None:
        cmd += ['--cut_ga']
    else:
        cmd += ['-E', str(evalue)]
    run(cmd + [rfam_db, fasta], verbose)
    hits = []
    for line in open(tblout):
        if line.startswith('#') or not line.strip():
            continue
        # target name, accession, query name, accession, mdl, mdl from, mdl to, seq from, seq to,
        # strand, trunc, pass, gc, bias, score, E-value, inc, description of target
        cols = line.split()
        hits.append({'family': cols[0], 'accession': cols[1], 'name': cols[2],
                     'seq_from': int(cols[7]), 'seq_to': int(cols[8]), 'strand': cols[9],
                     'score': float(cols[14]), 'evalue': float(cols[15])})
    if verbose:
        for h in hits:
            print(' hit: %(name)s %(family)s %(accession)s %(seq_from)s-%(seq_to)s %(strand)s score: %(score)s E: %(evalue)s' % h)
    return hits


def get_best_scores(hits):
    """:returns: {seq name: {family: best score}}, only hits on the + strand"""
    scores = {}
    for h in hits:
        if h['strand'] != '+':
            continue
        fam = scores.setdefault(h['name'], {})
        fam[h['family']] = max(fam.get(h['family'], h['score']), h['score'])
    return scores


def find_common_family(scores, name1, name2):
    """Get the family hit by both sequences with the highest sum of scores, or None."""
    fams1 = scores.get(name1, {})
    fams2 = scores.get(name2, {})
    common = set(fams1) & set(fams2)
    if not common:
        return None
    return max(common, key=lambda f: fams1[f] + fams2[f])


def cmalign(family, seqs, rfam_db, workdir, output_fn, verbose=False):
    """Fetch the family CM from Rfam and align seqs to it, save the alignment (Stockholm) to output_fn.

    Whole sequences are aligned; residues outside of the family (e.g. flanks)
    end up in insert columns ('.' in RF) and are not used."""
    cm = os.path.join(workdir, family + '.cm')
    fasta = os.path.join(workdir, family + '.fa')
    run(['cmfetch', '-o', cm, rfam_db, family], verbose)
    write_fasta(seqs, fasta)
    run(['cmalign', '-o', output_fn, cm, fasta], verbose)
    return output_fn


def get_seq_names(fns):
    """Get unique names, safe for Stockholm/fasta, for files (basename without extension)."""
    names = []
    for fn in fns:
        name = re.sub(r'[^\w.\-]', '_', os.path.splitext(os.path.basename(fn))[0]) or 'seq'
        new = name
        i = 2
        while new in names:
            new = '%s_%s' % (name, i)
            i += 1
        names.append(new)
    return names
