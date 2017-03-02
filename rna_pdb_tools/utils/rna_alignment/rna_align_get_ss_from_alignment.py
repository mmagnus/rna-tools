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
    """Get parser of arguments"""
    #try:
    #parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    #except:
    parser = argparse.ArgumentParser()#description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--file', help="subsection of an alignment", required=True)
    return parser

def rfam_ss_notat_to_dot_bracket_notat(c):
    """Take (c)haracter and standardize ss"""
    if c in [',', '_']:
        return '.'
    if c == '<':
        return '('
    if c == '>':
        return ')'
    return c

def clean_seq(seq,ss):
    #print seq, ss
    # new
    nseq = ''
    nss = ''
    for i,j in zip(seq,ss):
        if i != '-': # gap
            #print i,j
            nseq += i
            nss += rfam_ss_notat_to_dot_bracket_notat(j)
    return nseq.strip(), nss.strip()

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
