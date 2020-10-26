#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is an improved version of the script that uses the Rfam MySQL database online interface (thanks @akaped for this idea) (so you need to be connected to the Internet, of course). Redirect the output to the file.

.. image:: ../pngs/species.png

.. warning :: This scripts needs mysql-connector-python-rf module to connect the Rfam MySQL server, so install it before using: ``pip install mysql-connector-python-rf``.

Example::

    $ rna_alignment_get_species.py RF00004.stockholm.stk
    # STOCKHOLM 1.0
    Sorex-araneus-(European-shrew)                     AUCGCU-UCU----CGGCC--UUU-U

Examples 2::

    $ rna_alignment_get_species.py u5_rfam_u5only.stk --verbose
    # STOCKHOLM 1.0
    #=GF WK U5_spliceosomal_RNA
    #=GF NC 39.90
    #=GF RT The spliceosomal snRNAs of Caenorhabditis elegans.
    #=GF TC 40.00
    #=GF RN [3]
    (...)
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

"""
from __future__ import print_function

from rna_tools.tools.rna_alignment.rna_alignment import RNAalignment
from rna_tools.Seq import RNASequence
import pandas as pd
import argparse
import urllib
import sys
import mysql.connector


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument('--id-width', type=int, default=50, help="define width of ids, trim species name when longer than this")
    parser.add_argument('--evo-mapping')
    parser.add_argument('--evo-mapping-default', action="store_true")
    parser.add_argument('--one', action="store_true")
    parser.add_argument('--osfn', help="cache file")
    parser.add_argument('--rfam', action="store_true")
    parser.add_argument("alignment", help="alignment")
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

def get_species(id, ocfn, rfam, verbose=False):
    """Get the name of a species and taxonomy using the `Rfam MySQL database`_ based on the accession database (id).

    Args:
        id: is an accession number taken from the alignment, e. g. ``AABX02000022.1``
        ocfn: is a cache file that can be used to cache the results
        verbose: be verbose

    Returns:
        os, oc: see below

    These files are OS (species name) and OC (taxonomy)::

        OS   Leishmania tarentolae
        OC   Eukaryota; Euglenozoa; Kinetoplastida; Trypanosomatidae; Leishmaniinae;
        OC   Leishmania; lizard Leishmania.

    now this is done with MySQL (thanks akaped for this idea)

    .. _Rfam MySQL database: https://docs.rfam.org/en/latest/database.html

    """
    # clean for AABX02000022.1/363025-363047 -> AABX02000022.1
    id = clean_id(id)

    # cache
    if ocfn:
        try:
            df = pd.read_csv(ocfn, index_col=0)
        except:
            df = pd.DataFrame()
        if id in df.index:
            os = df[df.index == id]['os'].item()
            return os, ''  #


    if not rfam:
        from Bio import Entrez, SeqIO
        #provide your own mail here
        Entrez.email = "A.N.Other@example.com"
        try:
            handle = Entrez.esearch(db="nuccore", term=id, retmax='200')
            record = Entrez.read(handle)
        except:
            return id, 'Too many requsts'

        ids = record['IdList']
        for i in ids:
            try:
                handle = Entrez.efetch(db="nucleotide", id=i, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
            except:
                return id, 'Too many requsts'
            os = record.annotations["organism"]
            return os.replace(' ', '-'), ''
        return id, 'failed'

    else:
        mycursor = mydb.cursor()
        cmd = 'select t1.species, t1.tax_string from taxonomy t1 inner join rfamseq t2 on t1.ncbi_id = t2.ncbi_id where t2.rfamseq_acc like "' + id + '%" limit 1;'
        if verbose: print(cmd)
        mycursor.execute(cmd)
        try:
            os, oc = mycursor.fetchone()
        except:
            os, oc = id, ''  # just return accession

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

# main
if __name__ == '__main__':
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

    parser = get_parser()
    args = parser.parse_args()

    a = args.alignment
    name_width = args.id_width

    if args.rfam:
           mydb = mysql.connector.connect(
           host="mysql-rfam-public.ebi.ac.uk",
           user="rfamro",
           password="",
           port="4497",
           database="Rfam",
        )

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

                ################################################################################
                os, oc = get_species(id, args.osfn, args.rfam, args.verbose)
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
                    print((group.replace(' ','-') + '|' + os.replace(' ', '-')[:name_width] + '|' + semi_clean_id(id) + ss).ljust(name_width), seq.strip())

                else:
                    # this is no group
                    # ech, this is super messy, sorry!
                    group = ''
                    print((os.strip().replace(' ', '-')[:name_width] + str('')).ljust(name_width), seq.strip())

                if args.one:
                    break

            elif '#=GC SS_cons' in l:
                ss = l.replace('#=GC SS_cons', '')
                print('#=GC SS_cons'.ljust(name_width), ss.strip())
            # OK, i'm not sure, if this should be RF_cons or RF
            elif '#=GC RF_cons' in l:
                ss = l.replace('#=GC SS_cons', '')
                ss = l.replace('#=GC RF_cons', '')
                print('#=GC RF_cons'.ljust(name_width), ss.strip())
            elif '#=GC RF' in l:
                ss = l.replace('#=GC RF', '')
                print('#=GC RF'.ljust(name_width), ss.strip())
            else:
                print(l.strip())
