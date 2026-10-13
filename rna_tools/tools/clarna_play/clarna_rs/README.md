# clarna_rs — ClaRNA in Rust

A port of the ClaRNA classifier (`rna_clarna_run.py -bp+stack`) to Rust: same classifier libraries
(`../ClaRNAlib/lib/*.json.gz`), same output format, ~20–40× faster.

    cargo build --release          # needs Rust (https://rustup.rs)
    target/release/clarna_rs a.pdb                 # = rna_clarna_run.py -ipdb a.pdb -bp+stack
    target/release/clarna_rs --batch a.pdb b.pdb   # each output preceded by "### <file>"
    target/release/clarna_rs --dump DIR a.pdb      # also DIR/a.outCR/.dist.tsv/.contacts.tsv (regression format)

Options: `-bps`, `-stack`, `-bp+stack` (default), `-thresh 0.6`, `--lib DIR`.

## Build on another machine

1. Install Rust (once; user-local, no root needed; macOS or Linux):

       curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
       . "$HOME/.cargo/env"        # or open a new shell

   (or `brew install rustup && rustup-init`, or your distro's `rustup` package)
2. Build in the repo:

       cd rna_tools/tools/clarna_play/clarna_rs
       cargo build --release       # first build downloads the crates (needs internet), ~1 min

3. The binary is `target/release/clarna_rs`. PyMOL4RNA picks it up from there automatically.

The path of the classifier libraries (`../ClaRNAlib/lib`) is compiled into the binary, so build it inside the
repo; if you copy the binary somewhere else, pass `--lib /path/to/ClaRNAlib/lib`. Python is not needed to run it.
For a machine without internet: run `cargo vendor` on a connected machine and copy the `vendor/` folder along.
PyMOL4RNA (`bb`, `bbb`, ...) uses `target/release/clarna_rs` when it is built, otherwise the Python ClaRNA.

| | Python ClaRNA | clarna_rs |
|---|---|---|
| 20 regression structures, one process | 29.8 s | 0.73 s |
| 6chr (621 nt) | 9.4 s | 0.41 s |
| u6 (106 nt) | 4.5 s | 0.27 s (mostly loading the library) |

## Validation

Against the golden outputs of the Python ClaRNA (`../regression`):

    cd ../regression
    ../clarna_rs/target/release/clarna_rs --dump /tmp/rs $(grep -v '^#' inputs.txt | sed 's|^|../../|') > /dev/null
    clarna_regression.py check --candidate /tmp/rs --score-tol 0.01

- 20/20 structures: the same close doublets, the same 8421 classifier hits (pair, type, description,
  reference, alt flag) and the same output lines (pairs, types, order).
- 10/20 match to 1e-4 in every score. In the others, 6 hits differ by up to 0.005 (printed scores in 23 of
  2698 output lines, by up to 0.003). All of them are stacking "fuzzy" clusters matched almost exactly
  (~1e-6 Å) in ClaRNA's own training structures (1evv, 1kxk, 1q9a): their full-score bound is ~0.001 Å,
  below the float32 noise of the Python ClaRNA itself (it superimposes in float32 LAPACK); clarna_rs
  computes in float64.

## Layout

| file | ported from |
|---|---|
| `src/pdb.rs` | Bio.PDB.PDBParser (models, chains, altlocs by occupancy, duplicates), `StrucFile.get_info_chains` |
| `src/geom.rs` | vector maths, Bio.SVDSuperimposer (Kabsch) |
| `src/residue.rs` | `distances.Residue`, `Doublet` (only the parameters the bp and stacking libraries use) |
| `src/program.rs` | `clarna.ContactProgram` and the `MyExpr*` classes |
| `src/clarna.rs` | `extract_close_doublets`, `_find_contacts_using_library`, `make_structure_dict`, `display_bbstack` |

ClaRNA quirks kept on purpose (see `../regression/readme.md`): doublet distance 0, both directions classified,
the residue-index check in `compute_close_doublets` that drops some pairs when the structure has non-A/C/G/U
residues, and model-local residue indices for multi-model files. Not ported: base-phosphate/base-ribose/other
libraries, RNAView/MC-Annotate/FR3D modes, the NEXT:* backbone atoms (unused by the bp and stacking libraries).
