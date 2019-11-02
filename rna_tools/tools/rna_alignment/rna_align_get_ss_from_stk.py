#!/usr/bin/env python

"""Process an alignment in the Stockholm format to get sequences and secondary structures:

Example::

    $ rna_align_get_ss_from_stk.py aligns/gmp_RF01786.stockholm_p75_pk.sto
    AAOX01000007.1/31274-31356
    AAGAAUAUAGAACACUGUGAUGAGCGGUUUUUAUUUGCACUUUAAACCGCUUGGAGUGACUAGUGCAGCCGGCCAAUGAUCUA
    .(((.(((.(..(((......((((((((...............))))))))...)))...............).))).))).
    CP000724.1/3560727-3560809
    AAAAAUGUAGAGCAAAUGAACUGCAGGUAUACAUGGACGCCUUAAACUGCAGGGAUGUAGUGGCGUAACCGACUAACAAUAUU
    ((.(.(((((.(((......((((((.(......[[...[....).))))))...)))...]...]..]...))).)).).))
    AACY023761929.1/1009-1091
    AUAAUUUGGUGGGCGUUGAUGUGCCCUUUGUAUCUGGUCGCUUGAGGGGUACGGAGCCAAUAGCGAAACCGCCGCCGUCAUAG
    .((...((((((((......((((((((.......[.[......))))))))...)))......]...]...)))))...)).

"""
from rna_tools.tools.rna_alignment.rna_alignment import clean_seq_and_ss, RNAalignment
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="subsection of an alignment")
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args()
    a = RNAalignment(args.file)
    for s in a:
        s.remove_gaps()
        print(s)
        print(s.seq)
        print(s.ss)
