#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import argparse
from Bio import SeqIO
import re

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--one",
                        action="store_true", help="break")

    parser.add_argument("--save-to-file-each-seq")
    parser.add_argument("--format", action="store_true")
    parser.add_argument("file", help="input fasta file")
    parser.add_argument("searchfor", help="")
    parser.add_argument("ofile", help="output fasta file")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    seqs = []
    for record in SeqIO.parse(args.file, "fasta"):
        if args.searchfor.lower() in record.description.lower():
            if args.format:
                desc = record.description
                # CM000663.2/12922554-12922660 Homo sapiens chr 1, GRCh38 reference primary assembly.
                parts = re.split('[/-]', desc)
                # ['CM000666.2', '88684954', '88684848 Homo sapiens chromosome 4, GRCh38 reference primary assembly.']
                start = parts[1]
                end = parts[2].split(' ')[0]
                chr = parts[2].split(' ')[4].replace(',', '')
                # Chr  X 116850945-116850839
                # desc = 'Chr ' + chr.rjust(2) + ' ' + start + '-' + end
                if int(start) < int(end):
                    # desc = 'Chr' + chr + '/' + start + '-' + end
                    desc = 'Chr' + chr + '_' + start + '__' + end
                else:
                    desc = 'Chr' + chr + '_' + end + '__' + start
                record.description = desc
                print(record.description)
                record.name = ''
                record.id = ''
                print(record)

            seqs.append(record)

            if args.save_to_file_each_seq:
                f = open(record.description.replace('/', '-').replace(',','-') + '.fa', 'w')
                f.write(str(record.seq))
                f.close()

        #print(record.seq)

        if args.one:
            break

    SeqIO.write(seqs, args.ofile, "fasta")
