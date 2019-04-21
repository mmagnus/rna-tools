#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_dca_mapping.py

This function is divided into three parts (and panels in the output):

1. map interactions on the gapped sequence,
2. map interactions on the sequence without gaps,
3. map interactions on the sequence that might have insertions (insertions are marked with lower-case letters).

Gseq and Seq are files in the Fasta format with a header, sequence, and secondary structure. Warning: the script is not tested (yet) with secondary structure with pseudoknots. It should work as well, but it is not tested yet. Secondary structures should be defined in the Vienna  ("dot-bracket") notation.

`--noshort` will remove interactions that the distance between them is less than 6. For example, an interaction like [1,3] will be removed.

`--noss` will remove all DCA interactions that come because of secondary structure::

   (..)
   [1,4] DCA interaction will be removed.

Q: your interactions start from 0 or 1?

For a DCA interaction file, ' ' is used for separator at the moment.

Required:  rna_tools.SecondaryStructure to parse secondary structures in the Vienna format (dot-bracket notation).

Example::

    [mm] tpp$ git:(master) ✗ rna_dca_mapping.py --gseq 2gdi_gapped_RmExtraPart_RedAligment.fa --seq 2gdi.fa --dca tpp.ec.txt_52_sselected.csv --noshort --offset 10 --noss
    interactions:
    [(5, 91), (6, 90), (7, 89), (8, 88), (9, 87), (10, 53), (11, 52), (12, 51), (13, 50), (15, 48), (16, 30), (17, 29), (18, 28), (18, 71), (19, 27), (20, 21), (20, 26), (21, 22), (21, 23), (21, 25), (22, 23), (22, 24), (23, 24), (24, 25), (26, 27), (36, 47), (37, 46), (37, 47), (38, 45), (39, 40), (39, 44), (40, 41), (40, 43), (41, 42), (42, 43), (43, 44), (47, 48), (54, 55), (55, 56), (56, 57), (60, 83), (61, 82), (63, 64), (63, 65), (65, 78), (66, 77), (67, 76), (68, 75), (69, 72), (69, 73), (69, 74), (72, 74)]
    pairs [[6, 92], [7, 91], [8, 90], [9, 89], [10, 88], [11, 54], [12, 53], [13, 52], [14, 51], [16, 38], [17, 31], [18, 30], [19, 29], [21, 28], [22, 27], [23, 26], [60, 85], [61, 84], [62, 83], [66, 79], [67, 78], [68, 77], [69, 76]]
    UNMAPPED SCORES ////////////////////////////////////////////////////////////////////////////////////////////////////////
    123456789112345678921234567893123456789412345678951234567896123456789712345678981234567899123456789
    -----GGACUCGGGGUGCCCUUCUUGAAGGCUGAGAAA----------UACCCGUAUCACCUGAUCUGGAUAAUGCCAGCGUAGGGAAGUUC------------
    -----(((((((((.((((.(((..))))))......)----------..)))).....(((...((((......))))...)))..)))))------------
                   x                                x                                                         [16, 49]
                      x                                                    x                                  [19, 72]
                       x       x                                                                              [20, 28]
                        x     x                                                                               [21, 27]
    MAPPED SCORES //////////////////////////////////////////////////////////////////////////////////////////////////////////
    Mapped Interactions:
    [[1, 77], [2, 76], [3, 75], [4, 74], [5, 73], [6, 39], [7, 38], [8, 37], [9, 36], [11, 34], [12, 26], [13, 25], [14, 24], [14, 57], [15, 23], [16, 17], [16, 22], [17, 18], [17, 19], [17, 21], [18, 19], [18, 20], [19, 20], [20, 21], [22, 23], [40, 41], [41, 42], [42, 43], [46, 69], [47, 68], [49, 50], [49, 51], [51, 64], [52, 63], [53, 62], [54, 61], [55, 58], [55, 59], [55, 60], [58, 60]]
    pairs [[6, 92], [7, 91], [8, 90], [9, 89], [10, 88], [11, 54], [12, 53], [13, 52], [14, 51], [16, 38], [17, 31], [18, 30], [19, 29], [21, 28], [22, 27], [23, 26], [60, 85], [61, 84], [62, 83], [66, 79], [67, 78], [68, 77], [69, 76]]
    123456789112345678921234567893123456789412345678951234567896123456789712345678981234567899123456789
    GGACUCGGGGUGCCCUUCUUGAAGGCUGAGAAAUACCCGUAUCACCUGAUCUGGAUAAUGCCAGCGUAGGGAAGUUC
    (((((((((.((((.(((..))))))......)..)))).....(((...((((......))))...)))..)))))
              x                      x                                                [11, 34]
                 x                                          x                         [14, 57]
                  x       x                                                           [15, 23]
                   x     x                                                            [16, 22]
    GGACUCGGGGUGCCCUUCUgcgUGAAGGCUGAGAAAUACCCGUAUCACCUGAUCUGGAUAAUGCCAGCGUAGGGAAGUUC
    (((((((((.((((.(((.....))))))......)..)))).....(((...((((......))))...)))..)))))
    FINAL MAPPING //////////////////////////////////////////////////////////////////////////////////
    GGACUCGGGGUGCCCUUCUgcgUGAAGGCUGAGAAAUACCCGUAUCACCUGAUCUGGAUAAUGCCAGCGUAGGGAAGUUC
    (((((((((.((((.(((.....))))))......)..)))).....(((...((((......))))...)))..)))))
              x                         x                                             [11, 37]
                 x                                             x                      [14, 60]
                  x          x                                                        [15, 26]
                   x        x                                                         [16, 25]
                    x     x                                                           [17, 23]
    [[21, 47], [24, 70], [25, 36], [26, 35], [27, 33]]
    draw_dists([[21, 47], [24, 70], [25, 36], [26, 35], [27, 33]])
    output file: tpp.ec.txt_52_sselected.csv_mapped.csv

"""

import pandas as pd
import sys
import numpy as np
import argparse

from rna_tools.SecondaryStructure import parse_vienna_to_pairs

def rna_dca_mapping(seqfn, gseqfn, file_interactions, noss, noshort, offset, mss, verbose):
    """This function is deviede into

    .. warning:: in the line that we load parameters, watch for sep argument that defines seperator of your file (line21)
    """
    v = verbose

    #
    # Load the data, open files and parse information
    #
    # final ss and seq, ungapped
    f = open(seqfn)
    header = f.readline().strip()  # get rid of header
    seq = f.readline().strip()
    ss = f.readline().strip()
    # gapped
    f = open(gseqfn)
    header = f.readline().strip()  # get rid of header
    gseq = f.readline().strip()
    gss = f.readline().strip()
    # DCA
    df = pd.read_csv(file_interactions,sep=" ")
    interactions = zip(df['i'].tolist(), df['j'].tolist())
    #
    # Show input
    #
    # [(38, 51), (7, 110), (37, 52) from the input
    interactions.sort()
    print 'interactions:\n', interactions
    #
    # Process unmapped scores on gaped sequence
    # I panel
    #
    pairs = parse_vienna_to_pairs(gss)[0]
    print 'pairs', pairs
    print 'UNMAPPED SCORES ' + '/' * len(gseq)
    print '123456789112345678921234567893123456789412345678951234567896123456789712345678981234567899123456789'
    print gseq
    print gss

    for i in interactions:
        # form 0 or from 1 ?! ## be careful here! scores starts from 0 or 1 !?
        ij = [0,0]
        ij[0] = i[0] + 1
        ij[1] = i[1] + 1
        i = ij

        if noshort:
            delta = i[1] - i[0]
            if delta < 6:
                continue

        # remove interactions from gaps
        if gseq[i[0] - 1] == '-' or gseq[i[1] - 1] == '-':
            continue

        if noss:
            if i in pairs:
                continue

        line_new = 'x'.rjust(i[0]) + 'x'.rjust(i[1] - i[0]) + str(i).rjust(len(gseq) - i[1] + 10)
        print line_new
    #
    # How this mapping works, kurwa?
    # II panel
    #
    print 'MAPPED SCORES //' + '/' * len(gseq)
    mapped_interactions = []
    for i in interactions:
        ij = i
        #ij = [0,0]
        #ij[0] = i[0] - 1
        #ij[1] = i[1] - 1
        if v: print 'ij:', ij
        if v: print gseq[ij[0]], gseq[ij[1]]

        # ok, here I test if this pair comes from gaps, . -> -
        if gseq[ij[0]] == '-' or gseq[ij[1]] == '-':
            if v: print 'Removed Interaction:', ij
        else:
            [a,b]=[ij[0] - gseq[:ij[0]].count('-') +1, ij[1] - gseq[:ij[1]].count('-') + 1,]
            if v: print i, '->', a, b
            mapped_interactions.append([a, b])

    print 'Mapped Interactions:\n', mapped_interactions
    if v:
        for i in mapped_interactions:
            print str(i)

    pairs = parse_vienna_to_pairs(gss)[0]
    print 'pairs', pairs
    print '123456789112345678921234567893123456789412345678951234567896123456789712345678981234567899123456789'
    print gseq.replace('-','')  # what is the gap character, - or . ?
    print gss.replace('-', '')
    mapped_interactions.sort()

    filtered_interactions = []
    pairs = parse_vienna_to_pairs(gss.replace('-', ''))[0]
    for i in mapped_interactions:
        if noshort:
            delta = i[1] - i[0]
            if delta < 6:
                continue

        # remove if there is a gap
        if seq[i[0] - 1] == '-' or seq[i[1] - 1] == '-':
            continue

        if noss:
            if i in pairs:
                continue

        line_new = 'x'.rjust(i[0]) + 'x'.rjust(i[1] - i[0]) + str(i).rjust(len(seq) - i[1] + 10)
        print line_new
    #
    # How to include a gap in the mapping?
    # III panel, the final
    #
    print seq
    print ss

    pairs = parse_vienna_to_pairs(ss.replace('-', ''))[0]

    def n_lower_chars(string):
        """ https://stackoverflow.com/questions/10953189/count-lower-case-characters-in-a-string """
        return sum(1 for c in string if c.islower())

    print 'FINAL MAPPING //' + '/' * len(seq)
    nmapped_interactions = []
    for i in mapped_interactions:
        ij = i
        if v: print 'ij:', ij
        if v: print seq[ij[0] - 1], seq[ij[1] - 1]
        # gap mapping
        a = ij[0] + n_lower_chars(seq[:ij[0]])
        b = ij[1] + n_lower_chars(seq[:ij[1]])
        if v: print i, '->', a,b
        nmapped_interactions.append([a,b])
    mapped_interactions = nmapped_interactions

    print(seq)
    print(ss)
    pairs = parse_vienna_to_pairs(ss)[0]
    nmapped_interactions = []
    for i in mapped_interactions:
        if noshort:
            delta = i[1] - i[0]
            if delta < 6:
                continue

        if seq[i[0] - 1] == '-' or seq[i[1] - 1] == '-':
            continue

        if noss:
            if i in pairs:
                continue

        line_new = 'x'.rjust(i[0]) + 'x'.rjust(i[1] - i[0]) + str(i).rjust(len(seq) - i[1] + 10)
        nmapped_interactions.append([i[0], i[1]])
        print line_new
        if mss:
            print(ss)
    mapped_interactions = nmapped_interactions

    if offset:
        # e.g. filter_interaction is [[18, 71]], if offset is 10 then it will give you [[28, 81]]
        nmapped_interactions = [[x[0] + offset, x[1] + offset] for x in mapped_interactions]
    mapped_interactions = nmapped_interactions

    print mapped_interactions
    print 'draw_dists(' + str(mapped_interactions) + ')'
    print 'output file:', file_interactions+"_mapped.csv"
    a = pd.DataFrame(list(mapped_interactions), columns=["i","j"])
    a.to_csv(file_interactions+"_mapped.csv",sep=" ")


def get_parser():
    parser = argparse.ArgumentParser()  # description=__doc__ , formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--seq', help='seq fn in Fasta format', required=True)
    parser.add_argument('--gseq', help='gapped sequence and secondary structure (like in the alignment used for DCA) in Fasta format', required=True)
    parser.add_argument('--dca', help='file with parsed interactions', required=True)
    parser.add_argument('--offset', help="offset", type=int)
    parser.add_argument('--noss', help='filter out ss from plot', action='store_true')
    parser.add_argument('--mss', help='ss every each line', action='store_true')
    parser.add_argument('--verbose', help='be verbose', action='store_true')
    parser.add_argument('--noshort', help='filter out short interactions, dist in seq < 6 nt', action='store_true')
    return parser

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    rna_dca_mapping(args.seq, args.gseq, args.dca, args.noss, args.noshort, args.offset, args.mss, args.verbose)
