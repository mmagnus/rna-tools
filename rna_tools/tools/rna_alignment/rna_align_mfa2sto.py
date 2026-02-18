#!/usr/bin/env python
"""Convert a multi-FASTA alignment to Stockholm format.

Usage::

    $ ./rna_align_mfa2sto.py alignment.fasta
    $ cat alignment.mfa | ./rna_align_mfa2sto.py - > alignment.sto

The Stockholm output keeps the sequence identifiers (first token of each
FASTA header) and emits a minimal, single block `# STOCKHOLM 1.0` entry.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help="multi-FASTA alignment file or '-' for stdin",
    )
    return parser.parse_args()


def sanitize_name(header: str, seen: set[str]) -> str:
    """Create a Stockholm-friendly identifier that is unique."""
    token = header.split()[0]
    token = re.sub(r"[^A-Za-z0-9_.-]", "_", token) or "seq"
    candidate = token
    suffix = 2
    while candidate in seen:
        candidate = f"{token}_{suffix}"
        suffix += 1
    seen.add(candidate)
    return candidate


def infer_output_path(input_path: str) -> str:
    """Map input file name to `<stem>.sto`."""
    p = Path(input_path)
    if p.suffix:
        return str(p.with_suffix('.sto'))
    return f"{input_path}.sto"


def confirm_overwrite(path: str) -> None:
    reply = input(f"Output file '{path}' exists. Overwrite? [y/N]: ").strip().lower()
    if reply not in {"y", "yes"}:
        raise SystemExit("Aborted to avoid overwriting existing file.")


def read_mfasta(handle: Iterable[str]) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: str | None = None
    seq_chunks: List[str] = []
    seen_names: set[str] = set()

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                name = sanitize_name(header, seen_names)
                records.append((name, "".join(seq_chunks)))
            header = line[1:].strip()
            seq_chunks = []
            continue
        if header is None:
            raise ValueError("Sequence data encountered before any FASTA header.")
        seq_chunks.append(line.replace(" ", ""))

    if header is not None:
        name = sanitize_name(header, seen_names)
        records.append((name, "".join(seq_chunks)))

    return records


def write_stockholm(records: Sequence[Tuple[str, str]], handle) -> None:
    if not records:
        raise ValueError("No sequences found in the provided multi-FASTA file.")

    handle.write("# STOCKHOLM 1.0\n")
    for name, seq in records:
        handle.write(f"{name}\t{seq}\n")
    handle.write("//\n")


def main() -> None:
    args = parse_args()

    if args.input == "-":
        source = sys.stdin
    else:
        source = open(args.input)

    try:
        records = read_mfasta(source)
    finally:
        if source is not sys.stdin:
            source.close()

    if args.input == "-":
        if not hasattr(args, "output") or not args.output:
            raise SystemExit("Reading from stdin requires specifying -o/--output or redirect stdout.")
        output_path = args.output
    else:
        output_path = getattr(args, "output", None) or infer_output_path(args.input)

    if output_path == "-":
        destination = sys.stdout
    else:
        if Path(output_path).exists():
            confirm_overwrite(output_path)
        destination = open(output_path, "w")

    try:
        write_stockholm(records, destination)
    finally:
        if destination is not sys.stdout:
            destination.close()


if __name__ == "__main__":
    main()
