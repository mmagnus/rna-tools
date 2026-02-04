#!/usr/bin/env python
"""Convert Stockholm RNA alignments to a two-column CSV file.
GCGCGGAAACAAUGAUGAAUGGGUUUAAAUUGGGCACUUGACUCAUUUUGAGUUAGUAGUGCAACCGACCGUGCU
"""
import argparse
import csv
import sys
from typing import Iterable, List, Tuple

from Bio import AlignIO
from pathlib import Path


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "stockholm",
        help="input Stockholm alignment file; if no extension is provided, .stk is assumed",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output CSV file (defaults to <input>.csv)",
    )
    parser.add_argument(
        "--strip-gaps",
        action="store_true",
        help="remove '-' characters; by default gaps are preserved",
    )
    parser.add_argument(
        "--strip-periods",
        action="store_true",
        help="also remove '.' characters, which are sometimes used for gaps",
    )
    return parser


def resolve_stockholm_path(arg: str) -> Tuple[Path, Path]:
    """Return the actual Stockholm file path and the raw argument path."""
    raw_path = Path(arg)
    candidates = [raw_path]
    if not raw_path.suffix:
        candidates.append(raw_path.with_suffix(".stk"))
        candidates.append(raw_path.with_suffix(".sto"))

    for candidate in candidates:
        if candidate.exists():
            return candidate, raw_path

    raise FileNotFoundError(f"Could not find Stockholm file for '{arg}'")


def stockholm_rows(
    path: Path, strip_gaps: bool = False, strip_periods: bool = False
) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    with open(path) as handle:
        for alignment in AlignIO.parse(handle, "stockholm"):
            for record in alignment:
                sequence = str(record.seq)
                if strip_gaps:
                    sequence = sequence.replace("-", "")
                if strip_periods:
                    sequence = sequence.replace(".", "")
                rows.append((record.id, sequence))
    return rows


def write_csv(rows: Iterable[Tuple[str, str]], out_handle) -> None:
    writer = csv.writer(out_handle)
    writer.writerow(["key", "sequence"])
    for key, sequence in rows:
        writer.writerow([key, sequence])


def main() -> int:
    parser = get_parser()
    args = parser.parse_args()

    stockholm_path, raw_arg_path = resolve_stockholm_path(args.stockholm)

    rows = stockholm_rows(
        stockholm_path, strip_gaps=args.strip_gaps, strip_periods=args.strip_periods
    )

    if args.output == "-":
        write_csv(rows, sys.stdout)
    else:
        if args.output:
            output_path = Path(args.output)
        else:
            parent = raw_arg_path.parent
            name = raw_arg_path.name
            if name.endswith((".stk", ".sto")):
                name = name.rsplit(".", 1)[0]
            if not name:
                name = stockholm_path.stem
            output_path = parent / f"{name}.csv"
        with open(output_path, "w", newline="") as out_handle:
            write_csv(rows, out_handle)

    return 0


if __name__ == "__main__":

    main()
