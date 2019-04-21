#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The output you simply get from the screen, save it it to a file.

Example::

    python rna_alignment_get_species.py ../test_data/u6-with-species.stk ../test_data/u2-with-species.stk

Questions: is ``#=GC RF_cons?`` or ``#=GC RF?``

"""
from __future__ import print_function

from rna_tools.utils.rna_alignment.rna_alignment import RNAalignment
from rna_tools.Seq import RNASequence

import argparse
import urllib2
import sys


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('--id-width', type=int, default=50)
    parser.add_argument('--evo-mapping')
    parser.add_argument('--evo-mapping-default', action="store_true")
    parser.add_argument('--one', action="store_true")
    parser.add_argument('--u5', action="store_true")
    parser.add_argument('--calc-energy', action="store_true")
    parser.add_argument("alignment")
    parser.add_argument("alignment2")
    return parser


def ungap(x):
    return x.replace('-', '')


def semi_clean_id(id):
    id = id.split('/')[0]
    return id


def clean_id(id):
    id = id.split('/')[0]
    id = id.split('.')[0]
    return id

def get_species(id, verbose=False):
    """
    OS   Leishmania tarentolae
    OC   Eukaryota; Euglenozoa; Kinetoplastida; Trypanosomatidae; Leishmaniinae;
    OC   Leishmania; lizard Leishmania.

    link:
    https://www.ebi.ac.uk/ena/data/view/AANU01000000&display=text&download=txt&filename=AANU01000000.txt
    """
    # clean for AABX02000022.1/363025-363047 -> AABX02000022.1
    id = clean_id(id)

    # download
    url = "https://www.ebi.ac.uk/ena/data/view/%s&display=text&download=txt&filename=tmp.txt" % id
    response = urllib2.urlopen(url)
    oc = ''
    os = ''
    for l in response:
        if l.startswith('OS'):
            os = l.replace('OS', '')
        if l.startswith('OC'):
            oc += l.replace('OC', '').strip()

    if not os:
        if verbose:
            print(id)
            print(url)
        return None, None

    return os.strip(), oc.strip()


class seq():
    def __init__(self, id, seq):
        id = id.split('/')[0]
        id = id.split('.')[0]
        self.id = id
        self.seq = seq


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    a = args.alignment
    name_width = args.id_width

    b = args.alignment2
    name_width = args.id_width

    alist = []
    alist_done = []
    for l in open(a):
          if not l.startswith('#') and not l.startswith('//'):
                try:
                    id, sequence = l.split()
                except:
                    print('Error in ', l)
                    sys.exit(1)

                if id not in alist_done:
                    alist.append(seq(id, sequence))
                    alist_done.append(id)

          if '#=GC SS_cons' in l:
              ss_cons_a = l.replace('#=GC SS_cons', '').strip()
          if '#=GC RF_cons' in l:
              rf_cons_a = l.replace('#=GC RF_cons', '').strip()

    blist = []
    blist_done = []
    for l in open(b):
          if not l.startswith('#') and not l.startswith('//'):
                try:
                    id, sequence = l.split()
                except:
                    print('Error in ', l)
                    sys.exit(1)

                if id not in blist_done:
                    blist.append(seq(id, sequence))
                    blist_done.append(id)

          if '#=GC SS_cons' in l:
              ss_cons_b = l.replace('#=GC SS_cons', '').strip()
          if '#=GC RF_cons' in l:
              rf_cons_b = l.replace('#=GC RF_cons', '').strip()

    sep = '!!!!!!!'
    out = '# STOCKHOLM 1.0\n'
    done = []
    for seqa in alist: # record a
        for seqb in blist:
            if seqa.id not in done:
                if seqa.id == seqb.id:
                    out += seqa.id.ljust(name_width) + seqa.seq + sep + seqb.seq + '%%%\n'
                    done.append(seqa.id)

    out += '#=GC SS_cons'.ljust(name_width) + ss_cons_a + sep + ss_cons_b + '\n'
    out += '#=GC RF_cons'.ljust(name_width) + rf_cons_a + sep + rf_cons_b + '\n'
    out += '//'
    print(out)
