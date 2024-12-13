#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from Bio import AlignIO
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("sequence", help="", default="") # nargs='+')
    parser.add_argument("stk", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    verbose = args.verbose
    
    # Load the alignment
    alignment = AlignIO.read(args.stk, "stockholm")
    alignment_name = args.stk.split('/')[-1]
    # Define the target sequence ID
    target_id = args.sequence

    
    # Find the target sequence
    target_seq = None
    for record in alignment:
        if record.id == target_id:
            target_seq = record.seq
            break

    if target_seq is None:
        print(f"Target sequence '{target_id}' not found in the alignment.")
    else:
        seq_len = 0
        cols_no_coverage = 0
        # Collect columns where the target sequence has no gaps
        filtered_columns = []
        for col_index in range(alignment.get_alignment_length()):
            if target_seq[col_index] != "-":
                column = [record.seq[col_index] for record in alignment]
                filtered_columns.append(column)
                seq_len += 1

        # Display filtered columns
        if verbose: print(f"Columns where the target sequence '{target_id}' has no gaps:")
        for col in filtered_columns:
            if verbose: print(col, len(col))
            # Count the number of '-'
            gap_count = col.count('-')
            if verbose: print(f"Number of '-': {gap_count}")
            ratio_per_col = gap_count/(len(col) - 1)
            if verbose: print(f'ratio of gaps: {ratio_per_col}') # -1 target_sequence
            if ratio_per_col > 0.5:
                if verbose: print('WARNING: more than 50% gaps in column', cols_no_coverage)
                cols_no_coverage +=1
            if verbose: print()

        ratio = round((1 - (cols_no_coverage / seq_len)) * 100)
        print(f"Query {target_id} length {seq_len} in the alignment {alignment_name} number of columns with low coverage (<50%) for the query sequence {cols_no_coverage}; coverge for the query: {ratio}%") 
