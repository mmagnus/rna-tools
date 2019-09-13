#!/usr/bin/env python

"""Calculate statistics of foldability on an alignment.

The tool uses ENTRANA [1] to calculate, what the authors called, foldability (column: "foldability") of a given sequence into a given secondary structure.

Next, MC-Fold [2] is executed to calculate free energy (column: "mcsym") on the sequence and the secondary structure obtained based on the alignment. The secondary structure is used as constraints.

The third used program is RNAfold from the Vienna package [3]. Also, in this case the secondary structure obtained with rna-tools from the RNA alignment is used as constraints, columns: "mfe" (minimum free energy), "mfess" (secondary structure for minimum free energy state), "cfe" (minimum free energy of centroid), "cfess" (secondary structure for centroid, "diversity" (ensemble diversity), "efe" (free energy of the thermodynamic ensemble), "efess" (secondary structure for the thermodynamic ensemble), "freq" (frequency of mfe structure in ensemble). RNAfold is also executed in with "--enforceConstraint" where the constraints are enforced. This run gives analogous values as the default RNAfold, to all RNAfold column "_enforce" is added.

The tool is able to calculate the distance Levenshtein (the difference between the two sequences)(column: "distance") from the target sequence and all sequence in the alignment to test if there is a bias in the accuracy towards the most similar sequences.

Another tool used from the Vienna package is RNAeval. The tool calculates free energy for a given sequence and secondary structure.

The accuracy is expressed as the median of core RMSD of 10% the lowest core RMSD models for the given sequences.

.. image :: ../../rna_tools/tools/rna_alignment/test_data/foldability/rp17_foldability.png

The correlations::

    accuracy             1.000000
    cfe                  0.653813
    foldability          0.622038
    mfe                  0.607340
    efe                  0.585077
    diversity            0.404350
    eval                 0.349499
    cfe_enforce          0.311744
    mfe_enforce          0.302973
    efe_enforce          0.280929
    distance             0.256870
    freq                 0.037037
    diversity_enforce    0.018429
    mcsym                0.017533
    freq_enforce        -0.037991
    length              -0.340809

The data:

We tested correlations between the above-mentioned statistics, and the highest correlation, 0.65 () was achieved to the centroid free energy calculated with RNAFold, which suggests that to some extent this metric could be used to pick sequence from the alignment to pick sequences that are more likely to fold.

However, this needs further investigation and the detailed analysis an all cases and more folded sequences.

1. Su C, Weir JD, Zhang F, Yan H, Wu T. ENTRNA: a framework to predict RNA foldability. BMC Bioinformatics. BioMed Central 2019
2. Parisien M, Major F. The MC-Fold and MC-Sym pipeline infers RNA structure from sequence data. Nature 2008;452:51-5
3. Lorenz R, Bernhart SH, Honer Zu Siederdissen C, Tafer H, Flamm C, Stadler PF, et al. ViennaRNA Package 2.0. Algorithms Mol Biol. BioMed Central; 2011;6:26-14.

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
