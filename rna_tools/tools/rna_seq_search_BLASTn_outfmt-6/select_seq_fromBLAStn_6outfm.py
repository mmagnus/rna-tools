#!/usr/bin/env python
"""
Parse to BLASTn tabular output format 6 to get sequences in fasta format from database.

BLASTn tabular output format 6 - Column headers:
qseqid **sseqid** pident length mismatch gapopen **qstart** **qend** sstart send evalue bitscore

author: A. Zyla - azyla

.. warning:: Tested only on fasta files! and requires Biopython (tested with v1.68)
"""


import logging
import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import warnings
# logger
logger = logging.getLogger()
handler = logging.StreamHandler()
logger.addHandler(handler)


def get_headers(blastn_out6):
    """Parse blastn_out6 to get list of headers and ranges of sequences
    :param blastn_out6: BLASTn tabular output format 6
    :return: list of headers
    """
    headers_list = []
    seq_start_list = []
    seq_end_list = []
    param_list = []

    for line in blastn_out6:
        new_line = line.split('\t')
        #print new_line
        headers_list.append(new_line[1])
        headers_list = list(headers_list)  # you can put a list(set(headers_list)) to remove duplicates
        #print headers_list
        seq_start_list.append(new_line[6])
        seq_start_list=list(seq_start_list)

        seq_end_list.append(new_line[7])
        seq_end_list=list(seq_end_list)
        #print seq_start_list


        param_list = list([headers_list, seq_start_list, seq_end_list])
        #print param_list
    return param_list


0
def get_records(sequences_file, headers,seq_start, seq_end):
    """Get a list of records from fasta sequence file according to list headers and range of sequence

    :param sequences_file: Biopython object file with sequences - SeqIO
    :param headers: List of headers
    :param seq_start: A starting position of matched sequence
    :param seq_end: An ending position of matched sequence
    :return: List of Biopython objects - SeqRecord
    """
    rec_list = []  #list of records
    for record in sequences_file:
        seqid = record.id
        for header,i,j in zip(headers, seq_start, seq_end):
            i = int(i) - 1 # start -1 because python is counting from 0 :)
            j = int(j) - 1 # end
                # print i,j
                # print record.seq
                # record.seq = record.seq[i:j]
            if header.strip() == seqid: #check if IDs are the same
                record.seq=record.seq[i:j] # check a range
                rec_list.append(record)

    return rec_list


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("Seqfn", help="sequences in the Fasta format. e.q nt_database")
    parser.add_argument("blastn_out6", help="blastn_out6")

    parser.add_argument("--outfn", help="output aln file (default: seqfn .fasta -> _out.fasta)")
    return parser


# main
if __name__ == '__main__':

    args = get_parser().parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    if not args.outfn:
        args.outfn = args.Seqfn.replace('.fasta', '_out.fasta')

    print("If your  database of sequence is big, patience you must have, my dear user :) ")

    Seqfn=args.Seqfn
    blastn_out6 = open(args.blastn_out6)
    sequences_file = SeqIO.parse(Seqfn, 'fasta') #load sequences file

    list_of_param = get_headers(blastn_out6)
    #print list_of_param[1]
    headers = list_of_param[0]
    #print headers

    seq_start = list_of_param[1]
    #print seq_start
    seq_end = list_of_param[2]


    sequences = get_records(sequences_file,headers, seq_start,seq_end) #
    if sequences == []:
        print("Warning: No sequences found!")
    else:
        #print sequences
        SeqIO.write(sequences, args.outfn, "fasta")
        print('DONE \nSaved to:', args.outfn)


