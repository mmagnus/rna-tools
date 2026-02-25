#!/usr/bin/env python3
"""Extract RNA sequence from a PDB file."""

import argparse
from pathlib import Path
from collections import OrderedDict


RESN_TO_BASE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "T": "U",
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "U",
    "DU": "U",
    "ADE": "A",
    "CYT": "C",
    "GUA": "G",
    "URI": "U",
    "THY": "U",
    "PSU": "U",
    "H2U": "U",
    "5MU": "U",
    "OMU": "U",
    "M2G": "G",
    "1MG": "G",
    "7MG": "G",
    "OMG": "G",
    "2MG": "G",
    "5MC": "C",
    "OMC": "C",
    "1MA": "A",
    "M2A": "A",
    "I": "A",
}


def resn_to_base(resn: str, keep_unknown: bool = False) -> str:
    r = resn.strip().upper()
    if r in RESN_TO_BASE:
        return RESN_TO_BASE[r]
    if len(r) == 1 and r in {"A", "C", "G", "U", "T"}:
        return "U" if r == "T" else r
    if keep_unknown:
        return "N"
    return ""


def parse_pdb_sequences(path: str, chain_filter=None, keep_unknown: bool = False):
    seq_by_chain = OrderedDict()
    seen = set()

    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                continue

            atom_name = line[12:16].strip()
            if atom_name not in {"C1'", "C1*", "P"}:
                continue

            altloc = line[16].strip()
            if altloc not in {"", "A", "1"}:
                continue

            resn = line[17:20].strip()
            chain = line[21].strip() or "_"
            if chain_filter and chain not in chain_filter:
                continue

            resi = line[22:26].strip()
            icode = line[26].strip()
            resid = (chain, resi, icode)
            if resid in seen:
                continue
            seen.add(resid)

            base = resn_to_base(resn, keep_unknown=keep_unknown)
            if not base:
                continue

            seq_by_chain.setdefault(chain, [])
            seq_by_chain[chain].append(base)

    return {ch: "".join(seq) for ch, seq in seq_by_chain.items()}


def main():
    parser = argparse.ArgumentParser(description="Read RNA sequence from a PDB file.")
    parser.add_argument("pdb", help="input PDB file")
    parser.add_argument(
        "-c",
        "--chain",
        action="append",
        dest="chains",
        help="chain to extract (can be given multiple times)",
    )
    parser.add_argument(
        "--single",
        action="store_true",
        help="print one concatenated sequence instead of per-chain output",
    )
    parser.add_argument(
        "--keep-unknown",
        action="store_true",
        help="keep unknown residues as N",
    )
    args = parser.parse_args()

    chain_filter = set(args.chains) if args.chains else None
    seqs = parse_pdb_sequences(args.pdb, chain_filter=chain_filter, keep_unknown=args.keep_unknown)

    if not seqs:
        raise SystemExit("No RNA residues found.")

    if args.single:
        print("".join(seqs[ch] for ch in seqs))
        return

    pdb_id = Path(args.pdb).stem.split("_")[0]
    for ch, seq in seqs.items():
        print(f">{pdb_id}_{ch}")
        print(seq)


if __name__ == "__main__":
    main()
