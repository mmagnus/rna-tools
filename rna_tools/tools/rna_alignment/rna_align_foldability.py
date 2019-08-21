#!/usr/bin/env python

"""Calculate statistics of foldability on an alignment.

Example::

    $ python rna_align_foldability.py test_data/gmp_ref.sto test_data/gmp_foldability.csv
                                           cfess_enforce  distance  diversity
    0  ((((((.(((......(((((((((..............)))))))...      1.00       3.96
    1  (((.(((((((..((......(((((.(((((.......))).......      0.69       5.56
    2                                                         0.73       3.84
    3                                                         0.69       5.92
    4  (((....((((.((((((...((((((((...............))...      0.73       7.49
    5  ((((.(((.(((..((..((((((((.................)))...      0.75       7.92
    6  .(((((((((.((((.....(((((((((...............))...      0.72       5.83
    7                                                         0.72       7.35
    8  ...((((((((.(((......((((((.((....(((...)))..)...      0.65       4.86

       diversity_enforce    efe  efe_enforce
    0               2.89 -14.77       -13.75
    1               3.70 -19.52       -18.25
    2               0.00 -15.41         0.00
    3               0.00 -13.55         0.00
    4               2.46  -8.58        -6.91
    5               6.37 -20.72       -20.08
    6               2.92 -11.87       -11.38
    7               0.00 -14.59         0.00
    8               3.83 -21.16       -20.64

                                                   efess
    0  ((((((..........(((((((((..............)))))))...
    1  {{..(((((((..((......(((((.(({((.......}}).......
    2  .......((((.(((..(...(((((((((..............))...
    3  .....((((((.(((.....(((((((((.......{{...,}..,...
    4  .......{(((.,{{{{,...((((((((...............))...
    5  {(((.(((.(((..((..((((((((,{{...,.....)}}..)))...
    6  .(((((((((.((((.....(((((((((...............))...
    7  .....{,.{{{.(((......((((((((((..........,.}))...
    8  ...{(((((((.(((......((((((.((....(((...)))..)...

                             ...                         length   mcsym
    0                        ...                           75.0  -39.73
    1                        ...                           85.0  -37.89
    2                        ...                           84.0  -35.40
    3                        ...                           86.0  -36.11
    4                        ...                           83.0  -37.37
    5                        ...                           80.0  -43.59
    6                        ...                           84.0  -42.95
    7                        ...                           84.0  -36.55
    8                        ...                           85.0  -43.58

                          mcsym comment
    0  energy best dynamics programming
    1                         BP energy
    2                         BP energy
    3                         BP energy
    4  energy best dynamics programming
    5  energy best dynamics programming
    6  energy best dynamics programming
    7                         BP energy
    8  energy best dynamics programming

                                                mcsym ss   mfe mfe_enforce
    0  ((((((.(((......((((((((................))))))... -13.9       -12.9
    1  (((.(((((((..((......(((((.((.................... -18.0       -17.3
    2  ....(..((((.(((......((((((((................)... -14.0         0.0
    3  .(...((((((.(((......((((((((.................... -12.0         0.0
    4  (((....((((.(((......((((((((...............))...  -7.2        -6.1
    5  ((((.(((.(((......((((((((.................)))... -18.6       -18.6
    6  .(((((((((.(((......((((((((.................)... -10.5       -10.5
    7  ...(.((.(((.(((......((((((((................)... -12.8         0.0
    8  ...((((((((.(((......((((((.(.................... -19.8       -19.8

                                                   mfess
    0  ((((((..........(((((((((..............)))))))...
    1  ((..(((((((..((......(((((.(((((.......))).......
    2  .......((((.(((......(((((((((..............))...
    3  .....((((((.(((.....(((((((((....................
    4  .....................((((((((...............))...
    5  ((((.(((.(((..((..(((((((((((.........)))..)))...
    6  .(((((((((.((((.....(((((((((...............))...
    7  .....((.(((.(((......((((((((((............)))...
    8  ...((((((((.(((......((((((.((....(((...)))..)...

                                           mfess_enforce
    0  ((((((.(((......(((((((((..............)))))))...
    1  (((.(((((((..((......(((((.(((((.......))).......
    2                                              error
    3                                              error
    4  (((....((((.((((((...((((((((...............))...
    5  ((((.(((.(((..((..(((((((((((.........)))..)))...
    6  .(((((((((.((((.....(((((((((...............))...
    7                                              error
    8  ...((((((((.(((......((((((.((....(((...)))..)...

                                                     seq
    0  GCGCGGAAACAAUGAUGAAUGGGUUUAAAUUGGGCACUUGACUCAU...
    1  CUGUCGAAGAGACGCGAUGAAUCCCGCCCUGUAAUUCGGGCACCUC...
    2  AAUCAAUAGGGAAGCAACGAAGCAUAGCCUUUAUAUGGACACUUGG...
    3  AAAUAUUAUAGAGAUGUUGAAGUAUAUUCUAUUAUUGGGCACCUUA...
    4  AUUUUAAGAGGAAAUUUUGAACUAUAUACUUAUUUGGGCACUUUGU...
    5  UGCAAUGGGUGUGAUGAAGUCCGGACAGUAAUGUGGGCACUUAGUC...
    6  AAUAUUUUAGAAACUGAGAAGUAUAUCUUAUUAUUGGGCAUCUGGA...
    7  AUAACGGCACGAAGCAAUGAAAUGUUCGAUGUAACCGGGCACCUAU...
    8  AAAUUAAGGGGAAGCGUUGAGCCGCUACCCAUAUGUGGUUCACUCG...

                                                      ss
    0  ((((((.(((......((((((((................))))))...
    1  (((.(((((((..((......(((((.((....................
    2  ....(..((((.(((......((((((((................)...
    3  .(...((((((.(((......((((((((....................
    4  (((....((((.(((......((((((((...............))...
    5  ((((.(((.(((......((((((((.................)))...
    6  .(((((((((.(((......((((((((.................)...
    7  ...(.((.(((.(((......((((((((................)...
    8  ...((((((((.(((......((((((.(....................

    [9 rows x 26 columns]

"""
from __future__ import print_function
from rna_tools.Seq import RNASequence, load_fasta_ss_into_RNAseqs
from rna_tools.tools.rna_alignment.rna_alignment import clean_seq_and_ss, RNAalignment
import argparse
import pandas as pd

def hr(*text):
    print()
    print(' '.join(text))
    print('------------------------------')

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="an alignment in the Stokholm format, the first seq will be used to calculate distance to (#TODO pick any seq)")
    parser.add_argument('output', help="csv pandas file")
    parser.add_argument("--all-stars", help="this takes usully super long", action="store_true")
    parser.add_argument("--dev",  action="store_true")
    parser.add_argument("--skip-mcfold",  action="store_true")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
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

        dist = s.get_distance_to(target_seq)
        s.remove_gaps()
        length = len(s)

        s2 = RNASequence(s.seq, s.ss)
        cst = s.ss
        if args.all_stars:
            cst = s.ss.replace('.', '*')
            explore=1
        else:
            explore=''
        if not args.skip_mcfold:
            ss = list(s2.predict_ss(method='mcfold', constraints=cst, explore=explore, verbose=args.verbose))
            ss[0] = str(ss[0]) # str first value
            if args.verbose: print(ss)
        else:
            ss = '','', ''
        foldability = str(s2.get_foldability())
        eval = str(s2.eval(verbose = args.verbose))
        mcsym = ss[0]

        # RNAfoldX
        mfe, mfess, efe, efess, cfe, cfess, freq, diversity = s2.predict_ss("RNAfoldX", constraints=cst, enforce_constraint=False, verbose=args.verbose)

        mfe_enforce, mfess_enforce, efe_enforce, efess_enforce, cfe_enforce, cfess_enforce, freq_enforce, diversity_enforce = s2.predict_ss("RNAfoldX", constraints=cst, enforce_constraint=True, verbose=args.verbose)
        print('=======================================', s.id)
        print(s.seq)
        print(s.ss)
        print('foldability:', foldability)

        hr('RNAeval:', eval)
        print()
        print('MC-Fold')
        print(ss[1], ss[0], ss[2])

        hr('RNAfold with cst')
        print(cst, 'cst')
        print(mfess, mfe)
        print(efess, efe)
        print(cfess, cfe)
        print(freq, diversity)

        hr('RNAfold with enforce cst')
        print(cst, 'cst')
        print(mfess_enforce, mfe_enforce)
        print(efess_enforce, efe_enforce)
        print(cfess_enforce, cfe_enforce)
        print(freq_enforce, diversity_enforce)

        print # for space
        df = df.append({
            'id': s.id,
            'seq': s.seq,
            'ss': s.ss,
            'distance': dist,
            'foldability': foldability,
            'eval' : eval,
            'mcsym': ss[0],
            'mcsym ss': ss[1],
            'mcsym comment': ss[2],
            'length': length,

            'mfe' : mfe,
            'mfess' : mfess,
            'efe' : efe,
            'efess' : efess,
            'cfe' : cfe,
            'cfess' : cfess,
            'freq' : freq,
            'diversity' : diversity,

            'mfe_enforce' : mfe_enforce,
            'mfess_enforce' : mfess_enforce,
            'efe_enforce' : efe_enforce,
            'efess_enforce' : efess_enforce,
            'cfe_enforce' : cfe_enforce,
            'cfess_enforce' : cfess_enforce,
            'freq_enforce' : freq_enforce,
            'diversity_enforce' : diversity_enforce,

            }, ignore_index=True)

        if args.dev:
            break

    print(df)
    df.to_csv(args.output)
