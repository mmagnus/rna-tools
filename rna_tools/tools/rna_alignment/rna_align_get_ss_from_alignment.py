#!/usr/bin/env python

"""Input as a file::

  >ade
  GCU-U-CAUAUAAUCCUAAUGAUAUGG-UUUGGGA-GUUUCUACCAAGAG-CC--UUAAA-CUCUU---GAUUAUG-AAGU-
  (((.(.((((,,,(((((((_______.))))))).,,,,,,,,(((((((__.._____))))))...),,)))).)))).

to get::

  >ade
  GCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGU
  ((((((((...(((((((.......)))))))........((((((.......))))))..))))))))

"""
from rna_tools.tools.rna_alignment.rna_alignment import clean_seq_and_ss
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="subsection of an alignment")
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args()
    f = open(args.file)
    header = f.readline().strip()
    seq = f.readline()
    ss = f.readline()
    nseq, nss = clean_seq_and_ss(seq, ss)
    print(header)
    print(nseq)
    print(nss)
