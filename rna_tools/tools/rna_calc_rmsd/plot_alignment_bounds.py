#!/usr/bin/env python
"""Plot how many models start at each residue position.

magnus with codex"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Iterable, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot how many models start at each alignment position.",
    )
    parser.add_argument("csv_path", type=Path, help="Input CSV with alignment columns.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Path for the output PNG (defaults to <csv>.alignment_plot.png).",
    )
    parser.add_argument(
        "--length",
        type=int,
        default=1500,
        help="Axis limit for sequence length (default: 1500).",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=1,
        help="Option1: group counts into bins of this size (default: 1, i.e., no binning).",
    )
    parser.add_argument(
        "--cumulative",
        action="store_true",
        help="Option2: plot the cumulative number of models starting up to each position.",
    )
    return parser.parse_args()


def parse_alignment_value(value: str | None) -> int | None:
    """Extract the numeric residue index from tokens such as 'A:42'."""
    if not value:
        return None
    token = value.split(":")[-1]
    try:
        return int(token)
    except ValueError:
        return None


def collect_points(csv_path: Path) -> List[int]:
    starts: List[int] = []
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            start = parse_alignment_value(
                row.get("target_alignment_start") or row.get("alignment_start")
            )
            if start is None:
                continue
            starts.append(start)
    if not starts:
        raise SystemExit("No alignment start values found in the CSV.")
    return starts


def summarize_counts(counter: Counter[int], length: int, bin_size: int) -> tuple[list[float], list[int]]:
    bin_size = max(1, bin_size)
    if bin_size == 1:
        positions = list(range(1, length + 1))
        counts = [counter.get(pos, 0) for pos in positions]
        return positions, counts
    positions: list[float] = []
    counts: list[int] = []
    for start in range(1, length + 1, bin_size):
        end = min(start + bin_size - 1, length)
        positions.append((start + end) / 2.0)
        counts.append(sum(counter.get(pos, 0) for pos in range(start, end + 1)))
    return positions, counts


def to_cumulative(values: list[int]) -> list[int]:
    running = 0
    cumulative: list[int] = []
    for value in values:
        running += value
        cumulative.append(running)
    return cumulative


def plot_points(
    starts: Iterable[int],
    output: Path,
    length: int,
    bin_size: int,
    cumulative: bool,
) -> None:
    counter = Counter(starts)
    positions, counts = summarize_counts(counter, length, bin_size)
    if cumulative:
        counts = to_cumulative(counts)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(positions, counts, linewidth=1.5)
    ax.set_xlabel("Residue index (nt)")
    if bin_size == 1:
        ylabel = "# models starting at position"
    else:
        ylabel = f"# models starting per {bin_size}-nt bin"
    if cumulative:
        ylabel += " (cumulative)"
    ax.set_ylabel(ylabel)
    ax.set_title("Distribution of target alignment start positions")
    ax.set_xlim(1, length)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    print(f"Saved plot to {output}")


def main() -> None:
    args = parse_args()
    csv_path: Path = args.csv_path
    if not csv_path.exists():
        raise SystemExit(f"Input CSV not found: {csv_path}")
    output = args.output or csv_path.with_suffix(".alignment_plot.png")
    starts = collect_points(csv_path)
    plot_points(starts, output, args.length, args.bin_size, args.cumulative)


if __name__ == "__main__":
    main()
