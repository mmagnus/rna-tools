#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The output you simply get from the screen, save it it to a file.

Example::

    rna_alignment_get_species.py RF00004.stockholm.stk
    # STOCKHOLM 1.0
    Sorex-araneus-(European-shrew)                     AUCGCU-UCU----CGGCC--UUU-U

Examples 2::

    [dhcp177-lan203] Desktop$ rna_alignment_get_species.py u5_rfam_u5only.stk --verbose
    # STOCKHOLM 1.0
    #=GF WK U5_spliceosomal_RNA
    #=GF NC 39.90
    #=GF RT The spliceosomal snRNAs of Caenorhabditis elegans.
    #=GF TC 40.00
    #=GF RN [3]
    #=GF RM 2339054
    #=GF RL Nucleic Acids Res 1990;18:2633-2642.
    #=GF AU Gardner PP; 0000-0002-7808-1213
    #=GF CC methylation.
    #=GF CB cmcalibrate --mpi CM
    #=GF DR GO; 0046540; U4/U6 x U5 tri-snRNP complex;
    #=GF ID U5
    #=GF SS Published; PMID:2339054; Griffiths-Jones SR
    #=GF RA Thomas J, Lea K, Zucker-Aprison E, Blumenthal T
    #=GF SQ 180
    #=GF SM cmsearch --cpu 4 --verbose --nohmmonly -E 1000 -Z 549862.597050 CM SEQDB
    #=GF DE U5 spliceosomal RNA
    #=GF AC RF00020
    #=GF SE Zwieb C, The uRNA database, PMID:9016512; PMID:18390578
    #=GF GA 40.00
    #=GF BM cmbuild -F CM SEED
    #=GF TP Gene; snRNA; splicing;
    Bos-taurus-(cattle)                                GAUC-GUAUAAAUCUUUCGCCUUUUACUAAAGA-UUUCCG----UGG-A--GA-G
    Sorex-araneus-(European-shrew)                     GAUC-GUAUAAAUCUUUCGCCUUUUACUAAAGA-UUUCCG----UGG-A--GA-G
    Ictidomys-tridecemlineatus-(thirteen-lined-ground- GAUC-GUAUAAAUCUUUCGCCUUUUACUAAAGA-UUUCCG----UGG-A--GA-G
    Monodelphis-domestica-(gray-short-tailed-opossum)  GAUC-GUAUAAAUCUUUCGCCUUUUACUAAAGA-UUUCCG----UGG-A--GA-G
    Oryctolagus-cuniculus-(rabbit)                     GAUC-GUAUAAAUCUUUCGCCUUUUACUAAAGA-UUUCCG----UGG-A--GA-G
    Cavia-porcellus-(domestic-guinea-pig)              GAUC-GUAUAAAUCUUUCGCCUUUUACUAAAGA-UUUCCG----UGG-A--GA-G
    Ochotona-princeps-(American-pika)                  GAUC-GUAUAAAUCUUUCGCCUUUUACUAAAGA-UUUCCG----UGG-A--GA-G

.. note::

  This code has way more code than the name of the script says. This is customized script based on
  some script that did way more.
"""
from __future__ import print_function

from rna_tools.tools.rna_alignment.rna_alignment import RNAalignment
from rna_tools.Seq import RNASequence
import pandas as pd
import argparse
import urllib2

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument('--id-width', type=int, default=50)
    parser.add_argument('--evo-mapping')
    parser.add_argument('--evo-mapping-default', action="store_true")
    parser.add_argument('--one', action="store_true")
    parser.add_argument('--u5', action="store_true")
    parser.add_argument('--calc-energy', action="store_true")
    parser.add_argument('--osfn')
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

def get_species(id, ocfn, verbose=False):
    """
    OS   Leishmania tarentolae
    OC   Eukaryota; Euglenozoa; Kinetoplastida; Trypanosomatidae; Leishmaniinae;
    OC   Leishmania; lizard Leishmania.

    link:
    https://www.ebi.ac.uk/ena/data/view/AANU01000000&display=text&download=txt&filename=AANU01000000.txt
    """
    # clean for AABX02000022.1/363025-363047 -> AABX02000022.1
    id = clean_id(id)

    if ocfn:
        try:
            df = pd.read_csv(ocfn, index_col=0)
        except:
            df = pd.DataFrame()
        if id in df.index:
            os = df[df.index == id]['os'].item()
            return os, ''  #
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

    if ocfn:
        os, oc = os.strip(), oc.strip()
        df.index.name = 'id'
        df = df.append(pd.DataFrame({'os' : [os]}, index=[id])) # 'oc': [oc]
        df.to_csv(ocfn) # , index=False)
    return os, oc


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

    os_done = []

    cc = 1
    for l in open(a):
        if l.strip():
            if not l.startswith('#') and not l.startswith('//'):
                if args.debug:
                    print(cc)
                    cc += 1
                try:
                    id, seq = l.split()
                except:
                    print('Error in ', l)
                    sys.exit(1)
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
                os, oc = get_species(id, args.osfn)
                if not os:
                    os = id
                os = os.replace('.', '_') # remove dots from here
                # check if os if it's there already
                c = 1
                while 1:
                    if os not in os_done:
                        os_done.append(os) # Tupaia-chinensis-(Chinese-tree-shrew).1.2.3.4.5.6. fuck!
                        break
                    if len(os.split('.')) == 2:
                        os = os.replace('.' + str(c - 1), '.' + str(c))
                    else:
                        os += '.' + str(c)
                    if args.debug:
                        pass
                    c += 1
                ################################################################################
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
                    # this is no group
                    # ech, this is super messy, sorry!
                    group = ''
                    print((os.strip().replace(' ', '-')[:name_width] + str(energy)).ljust(name_width), seq.strip())

                if args.one:
                    break
                #else:
                #    print(os.replace(' ', '-') + str(energy)[:name_width].ljust(name_width), seq.strip())

            elif '#=GC SS_cons' in l:
                ss = l.replace('#=GC SS_cons', '')
                print('#=GC SS_cons'.ljust(name_width), ss.strip())

            # OK, i'm not sure, if this should be RF_cons or RF
            elif '#=GC RF_cons' in l:
                ss = l.replace('#=GC SS_cons', '')
                print('#=GC RF_cons'.ljust(name_width), ss.strip())
            elif '#=GC RF' in l:
                ss = l.replace('#=GC RF', '')
                print('#=GC RF'.ljust(name_width), ss.strip())
            else:
                print(l.strip())
