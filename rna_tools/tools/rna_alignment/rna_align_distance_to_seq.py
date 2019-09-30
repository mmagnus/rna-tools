#!/usr/bin/env python

"""Calculate

"Process an alignment in the Stockholm format to get sequences and secondary structures:

Example::

    $ rna_align_distance_to_seq.py test_data/gmp_ref.sto test_data/gmp_ref_distance.csv

       distance                            id
    0      1.00                           gmp
    1      0.69    AE000513.1/1919839-1919923
    2      0.73      BA000004.3/387918-388001
    3      0.69  ABFD02000011.1/154500-154585
    4      0.73      AE015927.1/474745-474827
    5      0.75                AAWL01000006.1
    6      0.72                    AM180355.1
    7      0.72      CP001116.1/102374-102457
    8      0.65    AJ965256.1/1260708-1260792

                                                     seq
    0  -----GCGCGGAAAC-AAUGAUGAAU--GGG-UUUA-AAUUGGGC-...
    1  CUGUCGAAGAGACGC-GAUGAAUCCC--GCC-CUGUAAUUCGGGC-...
    2  AAUCAAUAGGGAAGC-AACGAAGCAU--AGC-CUUU-AUAUGGAC-...
    3  AAAUAUUAUAGAGAU-GUUGAAGUAU--AUU-CUAUUA-UUGGGC-...
    4  AUUUUAAGAGGAAAU-UUUGAACUAU--AUA-CUU--AUUUGGGC-...
    5  --UGCAA-UGGGUGU-GAUGAAGUCC--GGA-CAGUAAUGUGGGC-...
    6  AAUAUUU-UAGAAAC-UGAGAAGUAU--AUC-UUAUUA-UUGGGC-...
    7  AUAACGGCACGAAGC-AAUGAAAUGU--UCG-AUGU-AACCGGGC-...
    8  AAAUUAAGGGGAAGC-GUUGAGCCGC--UAC-CCAU-AUGUGGUUC...

                                                     ss:
    0  (((((((((((.(((.......((((..(((.(................
    1  (((((((((((.(((.......((((..(((.(................
    2  (((((((((((.(((.......((((..(((.(................
    3  (((((((((((.(((.......((((..(((.(................
    4  (((((((((((.(((.......((((..(((.(................
    5  (((((((((((.(((.......((((..(((.(................
    6  (((((((((((.(((.......((((..(((.(................
    7  (((((((((((.(((.......((((..(((.(................
    8  (((((((((((.(((.......((((..(((.(................

"""
from rna_tools.tools.rna_alignment.rna_alignment import clean_seq_and_ss, RNAalignment
import argparse
import pandas as pd

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="an alignment in the Stokholm format, the first seq will be used to calculate distance to (#TODO pick any seq)")
    parser.add_argument('output', help="csv pandas file")
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args()
    a = RNAalignment(args.file)
    target_seq = ''
    df = pd.DataFrame()
    for s in a:
        if 'XX' in s.seq:  # skip line with X
            continue
        if not target_seq:
            target_seq = s.seq
        df = df.append({
            'id': s.id,
            'seq': s.seq,
            'ss:': s.ss,
            'distance': s.get_distance_to(target_seq),
            'length': len(s.seq.replace('-', '')),
            }, ignore_index=True)
    print(df)
    df.to_csv(args.output)
