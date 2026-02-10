#!/usr/bin/env python
"""Plot how many models start at each residue position.

magnus with codex"""
from __future__ import annotations

import argparse
import csv
import math
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
        default=1532,
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
    parser.add_argument(
        "--coverage",
        action="store_true",
        help="Plot residue coverage instead of counts of start positions.",
    )
    parser.add_argument(
        "--coverage-rmsd",
        "--coverge-rmsd",
        dest="coverage_rmsd",
        action="store_true",
        help="Plot residue coverage and per-residue average RMSD profile.",
    )
    parser.add_argument(
        "--coverage-rmsd-ends",
        action="store_true",
        help="Like --coverage-rmsd but the RMSD at residue i averages crops ending at i.",
    )
    parser.add_argument(
        "--report-position",
        type=int,
        help="Print coverage/RMSD statistics for the specified residue index.",
    )
    parser.add_argument(
        "--zoom",
        nargs=2,
        type=int,
        metavar=("START", "END"),
        help="Zoom the coverage axis to the specified residue interval (e.g., --zoom 1400 1532).",
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


def normalize_zoom(zoom_args: list[int] | tuple[int, int] | None, length: int) -> tuple[int, int] | None:
    """Clamp the requested zoom window to [1, length] and return it if valid."""
    if not zoom_args:
        return None
    start, end = zoom_args
    if start > end:
        start, end = end, start
    start = max(1, start)
    end = min(length, end)
    if start >= end:
        return None
    return (start, end)


def collect_points(csv_path: Path) -> List[int]:
    """Read alignment start positions (if present) from the CSV."""
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


def collect_alignment_ranges(csv_path: Path) -> list[tuple[int, int]]:
    """Return (start, end) tuples for every row with alignment bounds."""
    ranges: list[tuple[int, int]] = []
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            start = parse_alignment_value(
                row.get("target_alignment_start") or row.get("alignment_start")
            )
            end = parse_alignment_value(
                row.get("target_alignment_end") or row.get("alignment_end")
            )
            if start is None or end is None:
                continue
            low = min(start, end)
            high = max(start, end)
            ranges.append((low, high))
    if not ranges:
        raise SystemExit("No alignment start/end values found in the CSV.")
    return ranges


def parse_rmsd(value: str | None) -> float | None:
    if not value:
        return None
    try:
        rmsd = float(value)
    except ValueError:
        return None
    if math.isnan(rmsd):
        return None
    return rmsd


def collect_ranges_with_rmsd(csv_path: Path) -> list[tuple[int, int, float | None]]:
    """Return (start, end, rmsd) tuples, allowing rmsd to be None when absent."""
    ranges: list[tuple[int, int, float | None]] = []
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            start = parse_alignment_value(
                row.get("target_alignment_start") or row.get("alignment_start")
            )
            end = parse_alignment_value(
                row.get("target_alignment_end") or row.get("alignment_end")
            )
            if start is None or end is None:
                continue
            low = min(start, end)
            high = max(start, end)
            rmsd = parse_rmsd(row.get("rmsd"))
            ranges.append((low, high, rmsd))
    if not ranges:
        raise SystemExit("No alignment start/end values found in the CSV.")
    return ranges


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


def compute_coverage(
    ranges: Iterable[tuple[int, int]], length: int
) -> list[int]:
    """Compute residue coverage counts via a difference-array sweep."""
    if length < 1:
        return []
    diff = [0] * (length + 3)
    for start, end in ranges:
        if end < 1 or start > length:
            continue
        truncated_start = max(1, start)
        truncated_end = min(length, end)
        if truncated_end < truncated_start:
            continue
        diff[truncated_start] += 1
        diff[truncated_end + 1] -= 1
    coverage: list[int] = []
    running = 0
    for pos in range(1, length + 1):
        running += diff[pos]
        coverage.append(max(0, running))
    return coverage


def compute_coverage_and_rmsd(
    ranges: Iterable[tuple[int, int, float | None]], length: int
) -> tuple[list[int], list[float]]:
    """Return coverage counts and per-residue average RMSD (covering crops)."""
    if length < 1:
        return [], []
    cov_diff = [0] * (length + 3)
    rmsd_sum_diff = [0.0] * (length + 3)
    rmsd_count_diff = [0] * (length + 3)
    for start, end, rmsd in ranges:
        if end < 1 or start > length:
            continue
        truncated_start = max(1, start)
        truncated_end = min(length, end)
        if truncated_end < truncated_start:
            continue
        cov_diff[truncated_start] += 1
        cov_diff[truncated_end + 1] -= 1
        if rmsd is not None:
            rmsd_sum_diff[truncated_start] += rmsd
            rmsd_sum_diff[truncated_end + 1] -= rmsd
            rmsd_count_diff[truncated_start] += 1
            rmsd_count_diff[truncated_end + 1] -= 1
    coverage: list[int] = []
    rmsd_profile: list[float] = []
    running_cov = 0
    running_sum = 0.0
    running_count = 0
    for pos in range(1, length + 1):
        running_cov += cov_diff[pos]
        running_sum += rmsd_sum_diff[pos]
        running_count += rmsd_count_diff[pos]
        coverage.append(max(0, running_cov))
        if running_count > 0:
            rmsd_profile.append(running_sum / running_count)
        else:
            rmsd_profile.append(math.nan)
    return coverage, rmsd_profile


def compute_rmsd_by_end(
    ranges: Iterable[tuple[int, int, float | None]], length: int
) -> list[float]:
    """Return per-residue averages using only crops whose alignment ends there."""
    if length < 1:
        return []
    sums = [0.0] * (length + 2)
    counts = [0] * (length + 2)
    for start, end, rmsd in ranges:
        if rmsd is None:
            continue
        if end < 1 and start < 1:
            continue
        truncated_end = max(1, min(length, max(start, end)))
        sums[truncated_end] += rmsd
        counts[truncated_end] += 1
    profile: list[float] = []
    for pos in range(1, length + 1):
        if counts[pos]:
            profile.append(sums[pos] / counts[pos])
        else:
            profile.append(math.nan)
    return profile


def report_position_stats(
    ranges: Iterable[tuple[int, int, float | None]], position: int
) -> None:
    """Print coverage count plus covering/ending RMSD summaries for a residue."""
    covering: list[tuple[int, int, float | None]] = []
    covering_rmsds: list[float] = []
    ending_rmsds: list[float] = []
    for start, end, rmsd in ranges:
        if position < start or position > end:
            continue
        covering.append((start, end, rmsd))
        if rmsd is not None:
            covering_rmsds.append(rmsd)
            if end == position:
                ending_rmsds.append(rmsd)
    print(f"Residue {position} stats")
    print(f"  coverage count: {len(covering)}")
    if covering_rmsds:
        avg_covering = sum(covering_rmsds) / len(covering_rmsds)
        print(f"  avg RMSD (covering): {avg_covering:.4f} Å")
    else:
        print("  avg RMSD (covering): n/a")
    if ending_rmsds:
        avg_endings = sum(ending_rmsds) / len(ending_rmsds)
        print(
            f"  avg RMSD (ending exactly at residue): {avg_endings:.4f} Å"
        )
    else:
        print("  avg RMSD (ending exactly at residue): n/a")
    if covering_rmsds:
        formatted = ", ".join(f"{value:.2f}" for value in covering_rmsds)
        print(f"  RMSD values (covering): {formatted}")


def plot_points(
    starts: Iterable[int],
    output: Path,
    length: int,
    bin_size: int,
    cumulative: bool,
) -> None:
    starts_list = list(starts)
    counter = Counter(starts_list)
    positions, counts = summarize_counts(counter, length, bin_size)
    if cumulative:
        counts = to_cumulative(counts)
    peak_idx, peak_count = max(enumerate(counts), key=lambda item: item[1])
    peak_pos = positions[peak_idx]
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
    ax.set_title(f"Distribution of alignment starts (n={len(starts_list)})")
    ax.set_xlim(1, length)
    ax.grid(True, alpha=0.3)
    peak_y = peak_count * 1.05 + 0.05 * max(1, peak_count)
    ax.text(
        peak_pos,
        peak_y,
        f"peak @ {peak_pos:.0f}: {peak_count}",
        fontsize=9,
        ha="center",
        va="bottom",
    )
    fig.subplots_adjust(left=0.08, right=0.95, bottom=0.15, top=0.9, wspace=0.28)
    fig.savefig(output, dpi=200)
    print(f"Saved plot to {output}")


def plot_coverage(
    coverage: list[int],
    output: Path,
    length: int,
    sample_count: int,
    zoom: tuple[int, int] | None = None,
) -> None:
    """Plot the coverage trace and optional annotations on a single axis."""
    positions = list(range(1, length + 1))
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(positions, coverage, linewidth=1.5)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Coverage (# models)")
    ax.set_title(f"Residue coverage from alignments (n={sample_count})")
    x_min, x_max = 1, length
    if zoom is not None:
        x_min, x_max = zoom
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    if coverage:
        end_idx = min(length - 1, len(coverage) - 1)
        annotate_x = length if zoom is None else min(zoom[1], length)
        end_y = coverage[end_idx]
        ax.annotate(
            f"end x={annotate_x}",
            xy=(annotate_x, end_y),
            xytext=(-8, 10),
            textcoords="offset points",
            ha="right",
            va="bottom",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="none"),
        )
    if length >= 1500 and (zoom is None or (zoom[0] <= 1500 <= zoom[1])):
        ax.axvline(1500, color="tab:purple", linestyle="--", linewidth=1.0)
        coverage_at_1500 = ""
        idx_1500 = 1500 - 1
        if 0 <= idx_1500 < len(coverage):
            coverage_at_1500 = f"\ncoverage={coverage[idx_1500]}"
        ax.text(
            1500,
            ax.get_ylim()[1],
            f"x=1500{coverage_at_1500}",
            color="tab:purple",
            fontsize=8,
            ha="right",
            va="top",
            rotation=90,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.7, ec="none"),
        )
    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.15, top=0.9)
    fig.savefig(output, dpi=200)
    print(f"Saved coverage plot to {output}")


def plot_coverage_with_rmsd(
    coverage: list[int],
    rmsd_profile: list[float],
    output: Path,
    length: int,
    sample_count: int,
    zoom: tuple[int, int] | None = None,
) -> None:
    """Plot coverage on the primary axis with an RMSD profile on the twin axis."""
    positions = list(range(1, length + 1))
    fig, ax1 = plt.subplots(figsize=(11.5, 5))
    ax1.plot(positions, coverage, color="tab:blue", linewidth=1.4, label="Coverage")
    ax1.set_xlabel("Residue index")
    ax1.set_ylabel("Coverage (# models)", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    x_min, x_max = 1, length
    if zoom is not None:
        x_min, x_max = zoom
    ax1.set_xlim(x_min, x_max)
    ax1.grid(True, alpha=0.3)
    if coverage:
        end_idx = min(length - 1, len(coverage) - 1)
        annotate_x = length if zoom is None else min(zoom[1], length)
        end_y = coverage[end_idx]
        ax1.annotate(
            f"end x={annotate_x}",
            xy=(annotate_x, end_y),
            xytext=(-8, 10),
            textcoords="offset points",
            ha="right",
            va="bottom",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="none"),
        )
    ax2 = ax1.twinx()
    ax2.plot(
        positions,
        rmsd_profile,
        color="tab:red",
        linewidth=1.2,
        label="Average RMSD",
    )
    ax2.set_ylabel("Average RMSD (Å)", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax2.set_xlim(1, length)
    # Explain how the RMSD curve is derived (average of covering crops) near residue≈800.
    annotate_x = min(max(1, 800), length)
    finite_rmsd = [value for value in rmsd_profile if not math.isnan(value)]
    annotate_y = max(finite_rmsd) if finite_rmsd else 0.0
    ax2.text(
        annotate_x,
        annotate_y,
        "Avg RMSD(pos) = sum RMSD of crops covering pos / #covering",
        fontsize=9,
        color="tab:red",
        ha="center",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="none"),
    )
    if length >= 1500 and (zoom is None or (zoom[0] <= 1500 <= zoom[1])):
        ax1.axvline(1500, color="tab:purple", linestyle="--", linewidth=1.0)
        idx_1500 = 1500 - 1
        coverage_at_1500 = (
            f"\ncoverage={coverage[idx_1500]}"
            if 0 <= idx_1500 < len(coverage)
            else ""
        )
        value_at_1500 = math.nan
        if 0 <= idx_1500 < len(rmsd_profile):
            value_at_1500 = rmsd_profile[idx_1500]
        label = f"x=1500{coverage_at_1500}"
        if not math.isnan(value_at_1500):
            label += f"\nRMSD@1500={value_at_1500:.2f} Å"
        ax1.text(
            1500,
            ax1.get_ylim()[1],
            label,
            color="tab:purple",
            fontsize=8,
            ha="right",
            va="top",
            rotation=90,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.7, ec="none"),
        )
    title = (
        f"Residue coverage & RMSD profile (n={sample_count}, finite RMSD residues shown)"
    )
    fig.suptitle(title)
    fig.subplots_adjust(left=0.09, right=0.96, bottom=0.12, top=0.88)
    fig.savefig(output, dpi=200)
    print(f"Saved coverage+RMSD plot to {output}")


def main() -> None:
    args = parse_args()
    csv_path: Path = args.csv_path
    if not csv_path.exists():
        raise SystemExit(f"Input CSV not found: {csv_path}")
    if args.coverage_rmsd and args.coverage_rmsd_ends:
        raise SystemExit("Choose only one of --coverage-rmsd or --coverage-rmsd-ends.")
    zoom_range = normalize_zoom(args.zoom, args.length)
    ranges_with_rmsd: list[tuple[int, int, float | None]] | None = None
    if args.coverage_rmsd or args.coverage_rmsd_ends or args.report_position:
        ranges_with_rmsd = collect_ranges_with_rmsd(csv_path)
        if args.report_position:
            report_position_stats(ranges_with_rmsd, args.report_position)
    if args.coverage_rmsd:
        coverage, rmsd_profile = compute_coverage_and_rmsd(
            ranges_with_rmsd, args.length
        )
        rmsd_output = args.output or csv_path.with_suffix(".coverage_rmsd_plot.png")
        plot_coverage_with_rmsd(
            coverage,
            rmsd_profile,
            rmsd_output,
            args.length,
            len(ranges_with_rmsd),
            zoom_range,
        )
        return
    if args.coverage_rmsd_ends:
        coverage = compute_coverage(
            [(start, end) for start, end, _ in ranges_with_rmsd], args.length
        )
        rmsd_profile = compute_rmsd_by_end(ranges_with_rmsd, args.length)
        rmsd_output = args.output or csv_path.with_suffix(".coverage_rmsd_ends_plot.png")
        plot_coverage_with_rmsd(
            coverage,
            rmsd_profile,
            rmsd_output,
            args.length,
            len(ranges_with_rmsd),
            zoom_range,
        )
        return
    if args.coverage:
        if ranges_with_rmsd is not None:
            ranges = [(start, end) for start, end, _ in ranges_with_rmsd]
        else:
            ranges = collect_alignment_ranges(csv_path)
        coverage = compute_coverage(ranges, args.length)
        coverage_output = args.output or csv_path.with_suffix(".coverage_plot.png")
        plot_coverage(
            coverage,
            coverage_output,
            args.length,
            len(ranges),
            zoom_range,
        )
        return
    output = args.output or csv_path.with_suffix(".alignment_plot.png")
    starts = collect_points(csv_path)
    plot_points(starts, output, args.length, args.bin_size, args.cumulative)


if __name__ == "__main__":
    main()
