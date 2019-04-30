#!/usr/bin/env python
"""random_assignment_of_nucleotides.py - Random assignment of nucleotides for non-typical characters in the sequence alignment (arg --alignfn or fasta file with sequneces (arg --seqfn)

::

    R = G A (purine)
    Y = U C (pyrimidine)
    K = G U (keto)
    M = A C (amino)
    S = G C (strong bonds)
    W = A U (weak bonds)
    B = G U C (all but A)
    D = G A U (all but C)
    H = A C U (all but G)
    V = G C A (all but T)
    N = A G C U (any)

author: A. Zyla - azyla

.. warning:: Tested only on fasta files! and requires Biopython (tested with v1.68)
"""

import logging
import argparse
import random
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import warnings

# logger
logger = logging.getLogger()
handler = logging.StreamHandler()
logger.addHandler(handler)

def get_align(alignfn):
    """Get seq from an alignment with gaps.
    Args:
       alignfn (str): a path to an alignment
    Usage::
        >>> get_align('test_data/aln1.fasta')
        SingleLetterAlphabet() alignment with 2 rows and 13 columns
        AGGGGGACAGNYU 1
        CYGA------CGG 2
obj1', description='obj1', dbxrefs=[]), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])
    Returns:
         alignment
    """
    alignment = AlignIO.read(alignfn, 'fasta')
    #alignment = SeqIO.index(alignfn, 'fasta')
    return alignment

def get_sequences(seqfn):
    """Get seq from an fasta file.
    Args:
       seqfn (str): a path to a fasta file
    Usage::
        >>> get_align('test_data/fasta.fasta')

    Returns:
         [SeqRecord(seq=Seq('GGGYYGCCNRW', SingleLetterAlphabet()), id='1', name='1', description='1', dbxrefs=[]), SeqRecord(seq=Seq('GGRGYYGCCUURWAA', SingleLetterAlphabet()), id='1', name='1', description='1', dbxrefs=[])]

    """
    #alignment = AlignIO.read(alignfn, 'fasta')
    alignment = []
    for seq_record in SeqIO.parse(seqfn, 'fasta'):
        alignment.append(seq_record)
        #print alignment
    return alignment

def write_align(align, outfn):
    """Write cleaned alignment with Biopython.
    Args:
        align (obj): a cleaned alignment
        outfn (str): a path to a new alignment file
    Returns:
        none: writes to a file in fasta format
    """
    io = AlignIO
    io.write(alignment, outfn, 'fasta')

    logger.info('Alignment written to %s' % outfn)


def write_seq(seqfn, outfn):
    """Write cleaned alignment with Biopython.
    Args:
        align (obj): a cleaned alignment
        outfn (str): a path to a new alignment file
    Returns:
        none: writes to a file in fasta format
    """
    io = SeqIO
    io.write(alignment, outfn, 'fasta')

    logger.info('Alignment written to %s' % outfn)


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--alignfn", help="alignment in the Fasta format")
    parser.add_argument("--seqfn", help="sequences in the Fasta format")
    parser.add_argument("--outfn", help="output aln file (default: alnfn .fasta -> _out.fasta)")
    return parser

# main
if __name__ == '__main__':
    args = get_parser().parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    if not args.outfn:
        if args.alignfn:
            args.outfn = args.alignfn.replace('.fasta', '_out.fasta')

        if args.seqfn:
            args.outfn = args.seqfn.replace('.fasta', '_out.fasta')


    if args.alignfn:
        alignment=get_align(args.alignfn)
        for record in alignment:
            #print(record.seq)
            s = ''
            nseq = ''

            # replace R
            RR = ['R']
            R = ['G', 'A']
            new_R = random.choice(R)

            # replace Y
            YY = ['Y']
            Y = ['U', 'C']
            new_Y = random.choice(Y)

            # replace K
            KK = ['K']
            K = ['G', 'U']
            new_K = random.choice(K)

            # replace M
            MM = ['M']
            M = ['A', 'C']
            new_M = random.choice(M)

            # replace S
            SS = ['S']
            S = ['G', 'C']
            new_S = random.choice(S)

            # replace W
            WW = ['W']
            W = ['A', 'U']
            new_W = random.choice(W)

            # replace B
            BB = ['B']
            B = ['G', 'U', 'C']
            new_B = random.choice(B)

            # replace D
            DD = ['D']
            D = ['G', 'A', 'U']
            new_D = random.choice(D)

            # replace H
            HH = ['H']
            H = ['A', 'C', 'U']
            new_H = random.choice(H)

            # replace V
            VV = ['V']
            V = ['G', 'C', 'A']
            new_V = random.choice(V)

            # replace N
            NN = ['N']
            N = ['G', 'A', 'C', 'U']
            new_N = random.choice(N)

            # replace M
            MM = ['M']
            # replace T
            TT = ['T']

            for i in record.seq:
                if str(i) in RR:
                    i = i.replace('R', new_R)
                    nseq += i
                elif str(i) in YY:
                    i = i.replace('Y', new_Y)
                    nseq += i
                elif str(i) in KK:
                    i = i.replace('K', new_K)
                    nseq += i
                elif str(i) in SS:
                    i = i.replace('S', new_S)
                    nseq += i
                elif str(i) in WW:
                    i = i.replace('W', new_W)
                    nseq += i
                elif str(i) in BB:
                    i = i.replace('B', new_B)
                    nseq += i
                elif str(i) in DD:
                    i = i.replace('D', new_D)
                    nseq += i
                elif str(i) in HH:
                    i = i.replace('H', new_H)
                    nseq += i
                elif str(i) in VV:
                    i = i.replace('V', new_V)
                    nseq += i
                elif str(i) in NN:
                    i = i.replace('N', new_N)
                    nseq += i
                elif str(i) in MM:
                    i = i.replace('M', '-')
                    nseq += i
                elif str(i) in TT:
                    i = i.replace('T', 'U')
                    nseq += i
                else:
                    nseq += i
                seq = nseq  # cleaned sequence
                #
            # print seq, "XXX"
            # print type(record.seq)
            record.seq = Seq(seq)

    if args.seqfn:
        alignment=get_sequences(args.seqfn)

        for seq_record in alignment:
            #print seq_record.seq
            s = ''
            nseq = ''

            #replace R
            RR = ['R']
            R = ['G', 'A']
            new_R = random.choice(R)


            #replace Y
            YY = ['Y']
            Y = ['U', 'C']
            new_Y = random.choice(Y)

            #replace K
            KK =['K']
            K = ['G','U']
            new_K = random.choice(K)

            #replace M
            MM = ['M']
            M = ['A', 'C']
            new_M = random.choice(M)

            #replace S
            SS = ['S']
            S = ['G', 'C']
            new_S = random.choice(S)

            #replace W
            WW = ['W']
            W = ['A', 'U']
            new_W = random.choice(W)

            #replace B
            BB = ['B']
            B = ['G', 'U', 'C']
            new_B = random.choice(B)

            #replace D
            DD = ['D']
            D = ['G', 'A', 'U']
            new_D = random.choice(D)

            #replace H
            HH = ['H']
            H = ['A', 'C', 'U']
            new_H = random.choice(H)

            #replace V
            VV = ['V']
            V = ['G', 'C', 'A']
            new_V = random.choice(V)

            #replace N
            NN = ['N']
            N = ['G', 'A', 'C', 'U']
            new_N = random.choice(N)

            # replace M
            MM = ['M']
            # replace T
            TT=['T']

            for i in seq_record.seq:

                if str(i) in RR:
                    i = i.replace('R', new_R)
                    nseq += i
                elif str(i) in YY:
                    i = i.replace('Y', new_Y)
                    nseq += i
                elif str(i) in KK:
                    i = i.replace('K', new_K)
                    nseq += i
                elif str(i) in SS:
                    i = i.replace('S', new_S)
                    nseq += i
                elif str(i) in WW:
                    i = i.replace('W', new_W)
                    nseq += i
                elif str(i) in BB:
                    i = i.replace('B', new_B)
                    nseq += i
                elif str(i) in DD:
                    i = i.replace('D', new_D)
                    nseq += i
                elif str(i) in HH:
                    i = i.replace('H', new_H)
                    nseq += i
                elif str(i) in VV:
                    i = i.replace('V', new_V)
                    nseq += i
                elif str(i) in NN:
                    i = i.replace('N', new_N)
                    nseq += i
                elif str(i) in MM:
                    i = i.replace('M', '')
                    nseq += i
                elif str(i) in TT:
                    i = i.replace('T', 'U')
                    nseq += i
                else:
                    nseq += i
                seq = nseq  # cleaned sequence
                #print nseq
                #
            # print seq, "XXX"
            # print type(record.seq)
            seq_record.seq = Seq(seq)


    #print type(record.seq)
    #print type(alignment)
    align = alignment

    if args.alignfn:
        write_align(align, args.outfn)

    if args.seqfn:
        write_seq(align, args.outfn)


    print('DONE \nSaved to:', args.outfn)
