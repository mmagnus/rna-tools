#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyMOL plugin: annotate base pairs with NAPAIR (https://github.com/kubasvoboda4/NAPAIR)

Install
-------

1. NAPAIR (Python >= 3.10), in its own virtual env, which this plugin uses if present::

    git clone https://github.com/kubasvoboda4/NAPAIR.git
    cd NAPAIR
    python3 -m venv .venv
    .venv/bin/pip install -r requirements.txt

   (without .venv, NAPAIR's requirements must be installed for `python3` on PATH)

2. This plugin - load it once in PyMOL::

    PyMOL> run /path/to/pymol_napair.py

   or for every session, add to ~/.pymolrc::

    run /path/to/pymol_napair.py

   (PyMOL4RNA.py loads it automatically if it sits next to it)

3. Where is NAPAIR - nothing to do if NAPAIR is cloned next to this file's folder
   or next to the rna-tools repo (e.g. ~/src/rna-tools and ~/src/NAPAIR), or if
   this file is put inside the NAPAIR folder. Otherwise set NAPAIR_PATH. PyMOL
   started from the Dock/Finder does not see variables exported in ~/.zshrc,
   so set it in ~/.pymolrc, before the `run` line::

    python
    import os; os.environ['NAPAIR_PATH'] = '/path/to/NAPAIR'
    python end

   NAPAIR is looked up in: NAPAIR_PATH, assign_in_structure.py on PATH, the
   places above; if not found, `napair` prints the places it tried.

Usage
-----

::

    PyMOL> napair                  # all
    PyMOL> napair chain A, show=1
    PyMOL> help napair
"""
import inspect
import os
import shutil
import subprocess
import tempfile

from pymol import cmd

# under PyMOL's `run` __file__ points into PyMOL itself; the code object keeps the real path
_NAPAIR_HERE = os.path.abspath(inspect.currentframe().f_code.co_filename)


def _napair_exe(cmdline):
    """Run cmdline in a shell, return (stdout, stderr)."""
    p = subprocess.run(cmdline, shell=True, capture_output=True, text=True)
    return p.stdout.strip(), p.stderr.strip()

def _find_napair():
    """Return (NAPAIR folder or None, list of places tried).

    Looks at: NAPAIR_PATH env variable, assign_in_structure.py on PATH, the folder of this file
    and its parent, a NAPAIR checkout next to the rna-tools repo this file belongs to.
    """
    cands = []
    if os.environ.get('NAPAIR_PATH'):
        cands.append(os.environ['NAPAIR_PATH'])
    exe_path = shutil.which('assign_in_structure.py')
    if exe_path:
        cands.append(os.path.dirname(os.path.realpath(exe_path)))
    here = os.path.dirname(_NAPAIR_HERE)
    cands += [here] + [os.path.join(os.path.dirname(here), n) for n in ('NAPAIR', 'napair')]
    # <repo>/rna_tools/tools/PyMOL4RNA/pymol_napair.py -> <repo>/../NAPAIR
    repo = os.path.dirname(os.path.dirname(os.path.dirname(here)))
    if os.path.isdir(os.path.join(repo, 'rna_tools')):
        cands += [os.path.join(os.path.dirname(repo), n) for n in ('NAPAIR', 'napair')]
    tried = list(dict.fromkeys(cands))
    for c in tried:
        if os.path.isfile(os.path.join(c, 'assign_in_structure.py')):
            return c, tried
    return None, tried

# PyMOL's own mmCIF export leaves label_asym_id/label_seq_id empty and NAPAIR then finds no pairs,
# so save as PDB and let gemmi (from NAPAIR's env) write a proper mmCIF
NAPAIR_PDB2CIF = """import gemmi, sys
s = gemmi.read_structure(sys.argv[1])
s.setup_entities()
for m in s:
    for ch in m:
        for i, r in enumerate(ch):
            r.label_seq = i + 1
s.make_mmcif_document().write_file(sys.argv[2])
"""


def _dotbracket_pk(n, pairs):
    """Dot-bracket for pairs (0-based i < j), crossing pairs get [], {}, <>."""
    brackets = ['()', '[]', '{}', '<>']
    s = ['.'] * n
    levels = []  # pairs placed on each bracket level
    for i, j in sorted(pairs, key=lambda p: (p[1] - p[0])):
        if s[i] != '.' or s[j] != '.':
            continue  # base already used (e.g. a triple), keep the first, shortest pair
        for lvl, placed in enumerate(levels + [[]]):
            if not any(a < i < b < j or i < a < j < b for a, b in placed):
                if lvl == len(levels):
                    levels.append(placed)
                placed.append((i, j))
                if lvl < len(brackets):
                    s[i], s[j] = brackets[lvl][0], brackets[lvl][1]
                break
    return ''.join(s)


def napair(selection='', show=0, folder='', canonical_only=0):
    """Annotate base pairs of the selection (e.g. a chain) with NAPAIR (https://github.com/kubasvoboda4/NAPAIR).

    Prints the pairs (Leontis-Westhof family + subclass), the sequence and the dot-bracket
    of cWW pairs (pseudoknots as [], {}; chains separated by &). The selection is copied
    into a new object, napair_1, napair_2, ... (one per run), and the pairs are drawn on it.

    Example::

        PyMOL> napair
        PyMOL> napair chain A
        PyMOL> napair 1xjr and chain A, show=1

    Args:
        selection (str): PyMOL selection, typically one chain (default: all)
        show (int): 1 - show paired residues as sticks and draw pairs as dashes between C1'
                    (napair_N_cWW in gray, napair_N_nc in red)
        folder (str): where to keep input/output files (default: a new tmp folder)
        canonical_only (int): 1 - print only cWW pairs

    NAPAIR is found by _find_napair(); its .venv python is used if present, else python3.
    """
    import csv
    import glob
    NAPAIR_PATH, tried = _find_napair()
    if not NAPAIR_PATH:
        print('NAPAIR not found (tried: %s), set NAPAIR_PATH to the NAPAIR folder' % (', '.join(tried) or 'nothing'))
        return
    py = os.path.join(NAPAIR_PATH, '.venv', 'bin', 'python')
    if not os.path.isfile(py):
        py = shutil.which('python3') or 'python3'
    # e.g. `napair chain A` should not pick up the chain A of napair_* objects from earlier runs
    selection = selection.strip() or 'all'
    if 'napair_' not in selection and any(n.startswith('napair_') for n in cmd.get_names('objects')):
        selection = '(%s) and not napair_*' % selection
    folder = folder or tempfile.mkdtemp(prefix='napair_')
    os.makedirs(folder, exist_ok=True)
    name = 'napair'
    pdb, cif = os.path.join(folder, name + '.pdb'), os.path.join(folder, name + '.cif')
    cmd.save(pdb, selection)
    conv = os.path.join(folder, 'pdb2cif.py')
    with open(conv, 'w') as f:
        f.write(NAPAIR_PDB2CIF)
    # clean env: PyMOL.app leaks its own python paths into children
    env = 'env -i HOME=' + os.path.expanduser('~') + ' PATH=/usr/bin:/bin '
    out, err = _napair_exe(env + py + ' ' + conv + ' ' + pdb + ' ' + cif)
    if not os.path.isfile(cif):
        print('pdb -> cif conversion failed:', err)
        return
    cmdline = (env + py + ' ' + os.path.join(NAPAIR_PATH, 'assign_in_structure.py') + ' ' + cif
               + ' --boundaries ' + os.path.join(NAPAIR_PATH, 'parameter-boundaries.csv')
               + ' --out ' + folder + ' -f csv')
    print(cmdline)
    out, err = _napair_exe(cmdline)
    csvs = [fn for fn in glob.glob(os.path.join(folder, '*.csv'))]
    if not csvs:
        print(out)
        print(err.splitlines()[-1] if err else 'NAPAIR: no base pairs found')
        return

    # residues of the selection in order: (chain, resi) -> (index, resn)
    residues = []
    cmd.iterate('(%s) and polymer.nucleic and name C1\'' % selection,
                'residues.append((chain, resi, resn))', space={'residues': residues})
    index = {(ch, resi): k for k, (ch, resi, resn) in enumerate(residues)}

    pairs, seen = [], set()
    with open(csvs[0]) as f:
        for r in csv.DictReader(f):
            a = (r['auth_asym_id1'], r['auth_seq_id1'] + r['pdbx_PDB_ins_code1'])
            b = (r['auth_asym_id2'], r['auth_seq_id2'] + r['pdbx_PDB_ins_code2'])
            # one line per pair, the 5' residue first
            if index.get(b, 1e9) < index.get(a, 1e9):
                a, b = b, a
                n1, n2 = r['auth_comp_id2'], r['auth_comp_id1']
                fam = r['family'][0] + r['family'][2] + r['family'][1]  # cWH -> cHW
            else:
                n1, n2, fam = r['auth_comp_id1'], r['auth_comp_id2'], r['family']
            key = (a, b, fam)
            if key in seen:
                continue
            seen.add(key)
            pairs.append((a, b, n1, n2, fam, r['full_subclass']))
    pairs.sort(key=lambda p: (index.get(p[0], 1e9), index.get(p[1], 1e9)))

    for a, b, n1, n2, fam, sub in pairs:
        if int(canonical_only) and fam != 'cWW':
            continue
        print('%-3s %5s   %-3s %5s   %-4s %3s-%-3s  %s' % (a[0], a[1], b[0], b[1], fam, n1, n2, sub))
    ncww = sum(1 for p in pairs if p[4] == 'cWW')
    print('pairs: %d (cWW: %d, non-canonical: %d)' % (len(pairs), ncww, len(pairs) - ncww))

    # sequence + dot-bracket of cWW pairs, chains separated by &
    if residues:
        cww = [(index[a], index[b]) for a, b, _, _, fam, _ in pairs
               if fam == 'cWW' and a in index and b in index]
        dbn = _dotbracket_pk(len(residues), cww)
        letters = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', 'GTP': 'G', 'ATP': 'A', 'CTP': 'C', 'UTP': 'U', 'PSU': 'U'}
        seq, db = '', ''
        for k, (ch, resi, resn) in enumerate(residues):
            if k and ch != residues[k - 1][0]:
                seq, db = seq + '&', db + '&'
            seq += letters.get(resn, resn[-1] if len(resn) else 'N')
            db += dbn[k]
        first, last = residues[0], residues[-1]
        print('%s%s .. %s%s' % (first[0], first[1], last[0], last[1]))
        print(seq)
        print(db)
    print('csv: ' + csvs[0])

    # copy the whole selection into a new object: napair_1, napair_2, ...
    if int(canonical_only):
        pairs = [p for p in pairs if p[4] == 'cWW']
    k = 1
    while 'napair_%d' % k in cmd.get_names('all'):
        k += 1
    obj = 'napair_%d' % k

    def _res(r):
        return '(chain %s and resi \\%s)' % (r[0], r[1])
    cmd.create(obj, selection)
    print('object: ' + obj)

    if int(show):
        for a, b, _, _, fam, _ in pairs:
            cmd.show('sticks', '%s and (%s or %s)' % (obj, _res(a), _res(b)))
            cmd.distance(obj + ('_cWW' if fam == 'cWW' else '_nc'),
                         '%s and %s and name C1\'' % (obj, _res(a)), '%s and %s and name C1\'' % (obj, _res(b)))
        for o, color in ((obj + '_cWW', 'gray70'), (obj + '_nc', 'red')):
            if o in cmd.get_names('objects'):
                cmd.hide('labels', o)
                cmd.color(color, o)

cmd.extend('napair', napair)
print('napair - base pairs of a selection (e.g. chain A) with NAPAIR; napair chain A, show=1')

