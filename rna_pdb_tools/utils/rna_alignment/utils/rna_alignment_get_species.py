#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The output you simply get from the screen, save it it to a file.
"""
from __future__ import print_function

from rna_pdb_tools.utils.rna_alignment.rna_alignment import RNAalignment
from rna_pdb_tools.Seq import RNASequence

import argparse
import urllib2

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


# some default simple mapping
mapping =  [['Metazoa', 'Metazoa'],
           ['Viridiplantae', 'Viridiplantae'],
           ['Saccharomycotina', '(Fungi) Saccharomycotina'],
           ['Euglenozoa', 'Euglenozoa'],
           ['Alveolata', 'Alveolata'],
           ['Amoebozoa', 'Amoebozoa'],
           ['Cryptophyta', 'Cryptophyta'],
           ['Parabasalia', 'Parabasalia'],
           ['Rhizaria', 'Rhizaria'],
           ['Stramenopiles', 'Stramenopiles'],
           ['metagenomes', 'metagenomes'],
           ]


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    a = args.alignment
    name_width = args.id_width

    if args.evo_mapping:
        mapping = eval(open(args.evo_mapping).read().replace('\n', ''))
        if args.verbose:
            print(mapping)

    for l in open(a):
        if l:
            if not l.startswith('#') and not l.startswith('//'):
                id, seq = l.split()
                energy = ''
                ss = ''
                if False:
                    ################################################################################
                    a = ungap(seq[0:4])
                    b = ungap(seq[-3:])
                    ####################
                    # to the new alignment send only a + b, overwrite the original sequence
                    # print(seq)
                    seq = a + '-' + b
                    # print(seq)

                if args.u5:
                    a = ungap(seq[0:8])
                    b = ungap(seq[-9:])
                    seq = a + '-' + b

                if args.calc_energy:
                    # u5
                    seql = RNASequence(seq.replace('-','').replace('N', 'u'))
                    cst = "((((((((...........))))))))"
                    #AAAUCUUUCGCCUUUUACUAAAGA-UUU
                    #((((((((...........))))).)))
                    energy, ss = seql.predict_ss(method="mcfold", constraints=cst, verbose=args.verbose)
                    # u6
                    #pass
                    ## # u6atac
                    ## a = 'g' + ungap(seq[0:5])
                    ## b = ungap(seq[-4:]) + 'c'
                    ## loop = "guaa"
                    ## seql = RNASequence(a + loop + b)
                    ## cst = '((((.((..))))))'
                    ## if args.verbose:
                    ##     print('seq %s' % seq)
                    ##     print(a + loop + b)
                    ##     print(cst)
                    ## energy, ss = seql.predict_ss(method="mcfold", constraints=cst, verbose=args.verbose)
                    ## if args.verbose: print(energy, ss)
                ################################################################################
                os, oc = get_species(id)
                if not os:
                    os = id
                if args.evo_mapping or args.evo_mapping_default:
                    for m in mapping: # m is ['Parabasalia', 'Parabasalia']
                        if oc:
                            if m[0] in oc:
                                group = m[1]
                        else:
                            group = '???'
                    #
                    #print((os.replace(' ', '-')[:name_width] + '[' + group.replace(' ','-') + ']' + semi_clean_id(id) + '|' + str(round(energy,2)) +
                    #       '|' + a + '-' + b + ss).ljust(name_width), seq.strip()) # clean_id(id)
                    if energy:
                        print((os.replace(' ', '-')[:name_width] + '[' + group.replace(' ','-') + ']' + semi_clean_id(id) + '|' + str(round(energy,2)) +
                               '|' + ss).ljust(name_width), seq.strip()) # clean_id(id)
                    else:
                        print((group.replace(' ','-') + '|' + os.replace(' ', '-')[:name_width] + '|' + semi_clean_id(id) + ss).ljust(name_width), seq.strip())

                else:
                    print((os.replace(' ', '-')[:name_width] + '[' + group.replace(' ','-') + ']' + str(energy)).ljust(name_width), seq.strip())

                if args.one:
                    break

                #else:
                #    print(os.replace(' ', '-') + str(energy)[:name_width].ljust(name_width), seq.strip())

            elif '#=GC SS_cons' in l:
                ss = l.replace('#=GC SS_cons', '')
                print('#=GC SS_cons'.ljust(name_width), ss.strip())
            else:
                print(l.strip())
