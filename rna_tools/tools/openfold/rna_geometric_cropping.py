#!/usr/bin/env python3
"""
rna_geometric_cropping.py

As described in:

On inputs to deep learning for RNA 3D structure prediction
M Szikszai, M Magnus, S Kadyan, E Rivas
bioRxiv, 2025.02. 14.638364

Geometric cropping for RNA using a weighted graph + Dijkstra.

Inputs:
- sequence: RNA sequence (string)
- ss: dot-bracket secondary structure, same length as sequence

Graph:
- residues are vertices
- sequential neighbors have edges weight 1
- base-paired residues have edges weight 0

Crop:

- choose a seed residue
- run Dijkstra from seed
- take crop_size residues with smallest distance (ties resolved deterministically)

This is useful as a "spatial/structure-aware" crop that respects base pairing.
"""

from __future__ import annotations

import argparse
import heapq
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


class DotBracketError(ValueError):
    pass


def parse_dotbracket_pairs(ss: str) -> Dict[int, int]:
    """
    Parse dot-bracket into a dict i->j for paired residues (0-based indices).
    Supports (), [], {}, <> as independent bracket types.

    Returns:
        pairs: dict mapping index to its paired index (both directions).
    """
    bracket_pairs = {
        "(": ")",
        "[": "]",
        "{": "}",
        "<": ">",
    }
    opening = set(bracket_pairs.keys())
    closing = {v: k for k, v in bracket_pairs.items()}

    stacks: Dict[str, List[int]] = {op: [] for op in opening}
    pairs: Dict[int, int] = {}

    for i, ch in enumerate(ss):
        if ch in opening:
            stacks[ch].append(i)
        elif ch in closing:
            op = closing[ch]
            if not stacks[op]:
                raise DotBracketError(f"Unbalanced dot-bracket: extra '{ch}' at position {i}")
            j = stacks[op].pop()
            pairs[i] = j
            pairs[j] = i
        elif ch == ".":
            continue
        else:
            raise DotBracketError(f"Unsupported dot-bracket character '{ch}' at position {i}")

    leftovers = {op: st for op, st in stacks.items() if st}
    if leftovers:
        first_op = next(iter(leftovers))
        raise DotBracketError(f"Unbalanced dot-bracket: missing closing for '{first_op}'")

    return pairs


def build_weighted_adj_list(L: int, pairs: Dict[int, int]) -> List[List[Tuple[int, int]]]:
    """
    Build adjacency list with weights.
    Each entry adj[u] contains (v, w) edges.

    Sequential edges: weight 1
    Base-pair edges: weight 0
    """
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(L)]

    # sequential neighbors
    for i in range(L - 1):
        adj[i].append((i + 1, 1))
        adj[i + 1].append((i, 1))

    # base pairs
    for i, j in pairs.items():
        if i < j:
            adj[i].append((j, 0))
            adj[j].append((i, 0))

    return adj


def dijkstra(adj: List[List[Tuple[int, int]]], source: int) -> List[int]:
    """
    Dijkstra on nonnegative weights.
    Returns dist list (int), with large number for unreachable.
    """
    L = len(adj)
    INF = 10**12
    dist = [INF] * L
    dist[source] = 0
    pq: List[Tuple[int, int]] = [(0, source)]

    while pq:
        d, u = heapq.heappop(pq)
        if d != dist[u]:
            continue
        for v, w in adj[u]:
            nd = d + w
            if nd < dist[v]:
                dist[v] = nd
                heapq.heappush(pq, (nd, v))

    return dist


def geometric_crop_indices(
    sequence: str,
    ss: str,
    crop_size: int,
    seed: Optional[int] = None,
    rng_seed: Optional[int] = None,
) -> List[int]:
    """
    Compute one crop.

    Strategy:
    - pick seed (random if None)
    - run Dijkstra
    - choose crop_size residues with smallest (distance, index) for deterministic ties
    - return sorted indices

    Args:
        seed: 0-based seed index, or None to random choose
        rng_seed: for reproducible random seed selection
    """
    if len(sequence) != len(ss):
        raise ValueError("sequence and ss must have the same length")
    L = len(sequence)
    if not (1 <= crop_size <= L):
        raise ValueError("crop_size must be between 1 and sequence length")

    pairs = parse_dotbracket_pairs(ss)
    adj = build_weighted_adj_list(L, pairs)

    if seed is None:
        r = random.Random(rng_seed)
        seed = r.randrange(L)
    if not (0 <= seed < L):
        raise ValueError("seed out of range")

    dist = dijkstra(adj, seed)

    # pick smallest distances, break ties by index for reproducibility
    ranked = sorted(((dist[i], i) for i in range(L)))
    chosen = [i for _, i in ranked[:crop_size]]
    chosen.sort()
    return chosen


def mask_sequence(sequence: str, indices: List[int], mask_char: str = ".") -> str:
    """Return a masked sequence showing only chosen indices, others replaced by mask_char."""
    keep = set(indices)
    return "".join(sequence[i] if i in keep else mask_char for i in range(len(sequence)))


def mask_ss(ss: str, indices: List[int]) -> str:
    """
    Mask dot-bracket.
    Residues outside crop become '.'. If one side of a pair is missing, both become '.'.
    """
    keep = set(indices)
    pairs = parse_dotbracket_pairs(ss)

    ss_list = list(ss)

    # first blank everything outside
    for i in range(len(ss_list)):
        if i not in keep:
            ss_list[i] = "."

    # then fix broken pairs inside crop
    for i, j in pairs.items():
        if i < j:
            if (i in keep) ^ (j in keep):  # one in, one out
                ss_list[i] = "."
                ss_list[j] = "."
            # if both in keep, keep original bracket symbols

    # any non-dot outside keep already blanked
    return "".join(ss_list)


@dataclass
class CropResult:
    seed: int
    indices: List[int]
    seq_masked: str
    ss_masked: str


def generate_crops(
    sequence: str,
    ss: str,
    crop_size: int,
    n_crops: int = 5,
    rng_seed: int = 0,
) -> List[CropResult]:
    """
    Generate multiple crops by choosing different random seeds.
    """
    r = random.Random(rng_seed)
    results: List[CropResult] = []
    L = len(sequence)

    for _ in range(n_crops):
        seed = r.randrange(L)
        idx = geometric_crop_indices(sequence, ss, crop_size, seed=seed)
        results.append(
            CropResult(
                seed=seed,
                indices=idx,
                seq_masked=mask_sequence(sequence, idx, mask_char="."),
                ss_masked=mask_ss(ss, idx),
            )
        )
    return results


def main() -> None:
    p = argparse.ArgumentParser(description="Geometric cropping for RNA using Dijkstra on a weighted graph.")
    p.add_argument("--seq", required=True, help="RNA sequence")
    p.add_argument("--ss", required=True, help="Dot-bracket secondary structure, same length as seq")
    p.add_argument("--crop-size", type=int, required=True, help="Crop size (number of residues)")
    p.add_argument("--seed", type=int, default=None, help="Seed residue index (0-based). If omitted, random.")
    p.add_argument("--rng-seed", type=int, default=0, help="RNG seed used if --seed is omitted")
    p.add_argument("--n-crops", type=int, default=1, help="Generate N crops by random seeds (ignores --seed if >1)")
    args = p.parse_args()

    if args.n_crops > 1:
        crops = generate_crops(args.seq, args.ss, args.crop_size, n_crops=args.n_crops, rng_seed=args.rng_seed)
        for k, c in enumerate(crops, 1):
            print(f"> crop {k} seed={c.seed} size={len(c.indices)}")
            print("indices:", ",".join(map(str, c.indices)))
            print("seq:", c.seq_masked)
            print("ss: ", c.ss_masked)
            print()
    else:
        idx = geometric_crop_indices(args.seq, args.ss, args.crop_size, seed=args.seed, rng_seed=args.rng_seed)
        print("indices:", ",".join(map(str, idx)))
        print("seq:", mask_sequence(args.seq, idx, mask_char="."))
        print("ss: ", mask_ss(args.ss, idx))


if __name__ == "__main__":
    main()
