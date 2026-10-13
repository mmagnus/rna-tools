#!/usr/bin/env python
"""ClaRNA regression set: golden outputs of the Python ClaRNA, to check any change or port against.

For each structure in inputs.txt three files are kept in golden/:

    <name>.outCR         the usual ClaRNA text output (-bp+stack, score >= 0.6), as rna_clarna_run.py prints it
    <name>.dist.tsv      stage 1, close doublets: res_i, res_j, n_type, distance
                         (distance is always 0 in ClaRNA; the file pins down which doublets get classified)
    <name>.contacts.tsv  stage 2, all classifier hits (also < 0.6):
                         res_i, res_j, n_type, desc, full_desc, score, reference, is_alt

Usage::

    clarna_regression.py make                  # (re)write golden/ with the current Python ClaRNA
    clarna_regression.py check                 # recompute with the current Python ClaRNA, compare to golden/
    clarna_regression.py check --candidate DIR # compare files written by another implementation (same names/format)
    clarna_regression.py check --candidate DIR --score-tol 0.01   # scores may differ by 0.01, the rest exactly

Lines in .tsv files are sorted, floats are written with 6 decimals and compared with --tol (default 1e-4).
"""
from __future__ import print_function

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(os.path.dirname(HERE))  # rna_tools/tools
GOLDEN = os.path.join(HERE, 'golden')

CLARNA_OPTS = 'bp+stack'
MIN_SCORE = 0.6


def read_inputs():
    fns = []
    with open(os.path.join(HERE, 'inputs.txt')) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                fns.append(line)
    return fns


def name_of(fn):
    """1ehz_rpr.pdb -> 1ehz_rpr; names must be unique within inputs.txt."""
    return os.path.splitext(os.path.basename(fn))[0]


def run_python_clarna(fn):
    """Run the Python ClaRNA in-process (libraries cached), return (outcr, dist_lines, contact_lines)."""
    from rna_tools.tools.clarna_play import rna_clarna_run as r
    if CLARNA_OPTS not in r._cached_see_clarna:
        r.run_clarna_direct(fn, CLARNA_OPTS, MIN_SCORE)  # loads and caches the libraries
    see = r._cached_see_clarna[CLARNA_OPTS]
    g = see.cutils.start_clarna(fn)
    # the same as run_clarna_direct(), but from the graph we already have
    see.pdbclasslist = {}
    for (u, v, d) in g.edges(data=True):
        if d['type'] == 'contact' and d['weight'] > MIN_SCORE:
            see.pdbclasslist = see.make_structure_dict(u, v, d, see.pdbclasslist)
    outcr = 'Classifier: Clarna\nchains:  ' + r.StrucFile(fn).get_info_chains() + '\n' + see.display_bbstack(see.pdbclasslist)
    dist, contacts = [], []
    for u, v, d in g.edges(data=True):
        if d['type'] == 'dist':
            dist.append('%s\t%s\t%s\t%.6f' % (u, v, d['n_type'], d['weight']))
        elif d['type'] == 'contact':
            contacts.append('%s\t%s\t%s\t%s\t%s\t%.6f\t%s\t%d' % (
                u, v, d['n_type'], d['desc'], d['full_desc'], d['weight'], d['reference'], int(d['is_alt'])))
    return outcr, sorted(dist), sorted(contacts)


def write(name, outcr, dist, contacts, folder):
    with open(os.path.join(folder, name + '.outCR'), 'w') as f:
        f.write(outcr)
    with open(os.path.join(folder, name + '.dist.tsv'), 'w') as f:
        f.write('\n'.join(dist) + '\n')
    with open(os.path.join(folder, name + '.contacts.tsv'), 'w') as f:
        f.write('\n'.join(contacts) + '\n')


def read(name, folder):
    def lines(ext):
        with open(os.path.join(folder, name + ext)) as f:
            return [l for l in f.read().split('\n') if l]
    with open(os.path.join(folder, name + '.outCR')) as f:
        outcr = f.read()
    return outcr, lines('.dist.tsv'), lines('.contacts.tsv')


def _keyed(lines, nkey, fcol):
    """tsv lines -> {key columns: (float column, rest)}"""
    d = {}
    for l in lines:
        c = l.split('\t')
        d[tuple(c[:nkey])] = (float(c[fcol]), tuple(c[fcol + 1:]))
    return d


def _outcr_lines(text):
    """outCR -> [(pair part, [(type, score), ...])] for the pair lines"""
    res = []
    for l in text.strip().split('\n')[2:]:
        t = l[42:].split()
        res.append((l[:42], list(zip(t[0::2], [float(x) for x in t[1::2]]))))
    return res


def compare(name, gold, cand, tol, score_tol=None):
    """Return a list of differences (empty if the candidate matches the golden output).

    score_tol: tolerance for classifier scores (contacts and outCR); everything else must match exactly.
    """
    if score_tol is None:
        score_tol = tol
    errs = []
    g_out, g_dist, g_con = gold
    c_out, c_dist, c_con = cand
    # stage 1: close doublets (key: res_i, res_j, n_type)
    gd, cd = _keyed(g_dist, 3, 3), _keyed(c_dist, 3, 3)
    for k in sorted(set(gd) - set(cd)):
        errs.append('dist: missing %s' % (k,))
    for k in sorted(set(cd) - set(gd)):
        errs.append('dist: extra %s' % (k,))
    for k in sorted(set(gd) & set(cd)):
        if abs(gd[k][0] - cd[k][0]) > tol:
            errs.append('dist: %s %.6f != %.6f' % (k, gd[k][0], cd[k][0]))
    # stage 2: classifier hits (key: res_i, res_j, n_type, desc, full_desc; then score, reference, is_alt)
    gc, cc = _keyed(g_con, 5, 5), _keyed(c_con, 5, 5)
    for k in sorted(set(gc) - set(cc)):
        errs.append('contact: missing %s' % (k,))
    for k in sorted(set(cc) - set(gc)):
        errs.append('contact: extra %s' % (k,))
    for k in sorted(set(gc) & set(cc)):
        if abs(gc[k][0] - cc[k][0]) > score_tol:
            errs.append('contact: %s score %.6f != %.6f' % (k, gc[k][0], cc[k][0]))
        if gc[k][1] != cc[k][1]:
            errs.append('contact: %s %s != %s' % (k, gc[k][1], cc[k][1]))
    # final output: header, pairs, their order and types must match; scores within score_tol
    if g_out.strip().split('\n') != c_out.strip().split('\n') and score_tol > tol:
        gl, cl = _outcr_lines(g_out), _outcr_lines(c_out)
        same_shape = (g_out.split('\n')[:2] == c_out.split('\n')[:2] and len(gl) == len(cl) and
                      all(a[0] == b[0] and [t for t, _ in a[1]] == [t for t, _ in b[1]] for a, b in zip(gl, cl)))
        if same_shape:
            for a, b in zip(gl, cl):
                for (t, x), (_, y) in zip(a[1], b[1]):
                    if abs(x - y) > score_tol + 0.00005:  # + printed rounding
                        errs.append('outCR: %s %s %.4f != %.4f' % (a[0].strip(), t, x, y))
            return errs
    if g_out.strip().split('\n') != c_out.strip().split('\n'):
        gl, cl = g_out.strip().split('\n'), c_out.strip().split('\n')
        for l in sorted(set(gl) - set(cl)):
            errs.append('outCR: missing   %s' % l)
        for l in sorted(set(cl) - set(gl)):
            errs.append('outCR: extra     %s' % l)
        if not errs:
            errs.append('outCR: same lines, different order')
    return errs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mode', choices=['make', 'check'])
    p.add_argument('--candidate', help='folder with <name>.outCR/.dist.tsv/.contacts.tsv of another implementation')
    p.add_argument('--tol', type=float, default=1e-4)
    p.add_argument('--score-tol', type=float, help='tolerance for classifier scores (default: --tol); '
                   'pairs, types, references and line order must still match exactly')
    p.add_argument('--only', help='run only inputs whose name contains this string')
    args = p.parse_args()

    fns = [f for f in read_inputs() if not args.only or args.only in name_of(f)]
    names = [name_of(f) for f in fns]
    assert len(set(names)) == len(names), 'input names must be unique'

    # inputs that are not on this machine (e.g. large PDBs not in git) are skipped
    missing = [f for f in fns if not os.path.isfile(os.path.join(TOOLS, f))]
    for f in missing:
        print('%-55s SKIP (input not found: %s)' % (name_of(f), f))
    fns = [f for f in fns if f not in missing]

    if args.mode == 'make':
        if not os.path.isdir(GOLDEN):
            os.makedirs(GOLDEN)
        for fn in fns:
            t = time.time()
            outcr, dist, contacts = run_python_clarna(os.path.join(TOOLS, fn))
            write(name_of(fn), outcr, dist, contacts, GOLDEN)
            print('%-55s doublets %5d  hits %5d  %5.1fs' % (name_of(fn), len(dist), len(contacts), time.time() - t))
        return

    n_bad = 0
    for fn in fns:
        name = name_of(fn)
        t = time.time()
        if not os.path.isfile(os.path.join(GOLDEN, name + '.outCR')):
            print('%-55s SKIP (no golden output)' % name)
            continue
        gold = read(name, GOLDEN)
        if args.candidate:
            try:
                cand = read(name, args.candidate)
            except IOError as e:
                print('%-55s MISSING (%s)' % (name, e))
                n_bad += 1
                continue
        else:
            cand = run_python_clarna(os.path.join(TOOLS, fn))
        errs = compare(name, gold, cand, args.tol, args.score_tol)
        print('%-55s %s  %5.1fs' % (name, 'OK  ' if not errs else 'FAIL (%d)' % len(errs), time.time() - t))
        for e in errs[:10]:
            print('    ' + e)
        if len(errs) > 10:
            print('    ... %d more' % (len(errs) - 10))
        n_bad += bool(errs)
    n = len([f for f in fns if os.path.isfile(os.path.join(GOLDEN, name_of(f) + '.outCR'))])
    print('\n%d/%d structures match the golden output' % (n - n_bad, n))
    sys.exit(1 if n_bad else 0)


if __name__ == '__main__':
    main()
