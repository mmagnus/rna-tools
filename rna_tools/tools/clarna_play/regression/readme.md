# ClaRNA regression set

Golden outputs of the Python ClaRNA (`-bp+stack`, threshold 0.6) for the 20 structures in `inputs.txt`
(26–621 nt; cleaned and raw PDB entries, RNA-Puzzles models). `6chr*.pdb` (621 nt, 3.3 MB) and their golden files
are not in git; inputs or golden files that are missing are skipped, so 18 structures are checked elsewhere. Use it to check any change to ClaRNA,
and as the gate for a port (e.g. to Rust).

    clarna_regression.py make                    # rewrite golden/ (only when a change of results is intended)
    clarna_regression.py check                   # current Python ClaRNA vs golden/ (~30 s)
    clarna_regression.py check --candidate DIR   # another implementation's files vs golden/
    clarna_regression.py check --only 1ehz       # a subset
    clarna_regression.py check --candidate DIR --score-tol 0.01   # scores within 0.01, everything else exact

Run with the repo on the path, e.g. `PYTHONPATH=~/work/src/rna-tools ~/miniconda3/bin/python clarna_regression.py check`.

Per structure, `golden/` holds three files, so a port can be checked stage by stage:

| file | stage | content |
|---|---|---|
| `<name>.dist.tsv` | 1. close doublets | res_i, res_j, n_type, distance |
| `<name>.contacts.tsv` | 2. classifier | res_i, res_j, n_type, desc, full_desc, score, reference, is_alt — all hits, also < 0.6 |
| `<name>.outCR` | 3. output | the usual ClaRNA text output |

Scores are compared with tolerance 1e-4 (`--tol`); the text output must match line by line.
ClaRNA is deterministic: `check` reproduces 20/20 exactly.
The Rust port (`../clarna_rs`) passes 20/20 with `--score-tol 0.01` (see its README for why not all scores
match to 1e-4).

## Quirks a port must reproduce

- The doublet "distance" is always 0 (`close_doublets.append((i,j,0))`); the file only pins down which pairs are classified.
- Every close pair is classified in both directions, (i, j) and (j, i).
- Neighbours from the KD-tree are sorted (`utils.compute_close_doublets`). Before, scipy returned them in tree
  order, so the order of output lines with the same first residue depended on scipy internals; the golden
  files were made after this change (content unchanged, only line order).
- `compute_close_doublets` compares a residue index with the number of residues that have coordinates
  (`j >= len(residues_atoms_np)`): with waters or modified residues in the file, pairs with residues near the
  end of the chain are dropped (e.g. 1ehz: 358 doublets vs 748 in 1ehz_rpr).
- Residues other than A/C/G/U (after upper-casing) are kept in the residue list but skipped for doublets.
- Residues missing atoms print `ignoring ... (missing atoms)` to stdout (not part of `.outCR`).
- Alternative hits (`is_alt` = 1) are added when the best hit is uncertain (`?`) and the alternative scores > 0.7 × best.
