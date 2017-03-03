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
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--file', help="subsection of an alignment",  required=True)
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args() 

    f = open(args.file)
    header = f.readline().strip()
    seq = f.readline()
    ss = f.readline()
    nseq, nss = clean_seq(seq,ss)
    print header
    print nseq
    print nss
