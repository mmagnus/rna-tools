#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate energy (.cet) format::

    UGGC-CCCUGCGCAA-GGAUGACA
    (((..((((.....).))..))))
    (((..(((((***)).))..))))

Examples::

    $ rna_alignment_calc_energy.py --template alignments/u6-lower.cet alignments/u6-only-RemovedGapped.stk -v
          --loop-upper guaa --loop-lower guaa
          --loop-upper-cst '(..)' --loop-lower-cst '(..)'
    calc-energy2.py --template u6atac-template.txt u6atac_u6only.sto -v
   ./calc-energy2.py --template alignments/u6-lower.cet --one alignments/u6-lower-stem-only.sto

Takes cet files (calc-energy-templets)::

    $ rna_alignment_calc_energy.py --template test_data/u6-lower.cet --one test_data/u6-only.stk  -v # --loop-seq test_data/u6-only-loop-seq-u6-lower
    N/A% (0 of 182) |                                                                                                                            | Elapsed Time: 0:00:00 ETA:  --:--:--================================================================================
    AB010698.1/46467-46488
    (((..((((.....).))..))))
    UGGC-CCCUGCGCAA-GGAUGACA
    lower ------------------------------------
    UGG ugcgca ACA
    (((******)))
    UGGugcgcaACA
    (((((..))))) -10.64
    upper ------------------------------------
    UGGC-CCCUGCGCAA-GGAUGACA
    CCC ugcgca AGG
    CCCugcgcaAGG
    (((((..))))) -9.6
                           id  low_energy       low_seq        low_ss  up_energy        up_seq         up_ss
    0  AB010698.1/46467-46488      -10.64  UGGugcgcaACA  (((((..)))))       -9.6  CCCugcgcaAGG  (((((..)))))
    Done: u6-only-loop-seq-u6-lower

by parsing output from MC-Sym::

    domains have 5451 elements.
     10:47:16 up 141 days, 26 min,  0 users,  load average: 1.45, 1.30, 1.56
    Score: -999.000 GAACAUGGUUCUUGCCUUUUACCAGAACCAUCCGGGUGUUG
    Total number of MB structures with 3 stems: 16041
    (overlaps: 0, !energy: 335585)
    </pre><P><H2>Sorting the structures...
    <P></H2><pre></pre><H2><P><P><P>Filtered and Sorted solutions:<P><P><P></H2><pre>
    </pre><H2><P><P><P><a HREF="http://biwww2.informatik.uni-freiburg.de/Software/MARNA/index.html" target="_blank">MARNA</a>-formatted:<P><P><P></H2><pre>
    GAACAUGGUUCUUGCCUUUUACCAGAACCAUCCGGGUGUUG
    ((((((((((..))))))))))((((((((...)))))))) -33.20 ( -0.69)
    (((((((((....)))))))))((((((((...)))))))) -33.17 ( -0.69)
    ((((((((((((((((...))))))))))))......)))) -32.40 ( +0.00)

    Backtracking with 2 variables (stems):
    domains have 5451 elements.
     10:47:16 up 141 days, 26 min,  0 users,  load average: 1.45, 1.30, 1.56
    Score: -999.000 GAACAUGGUUCUUGCCUUUUACCAGAACCAUCCGGGUGUUG
    Total number of MB structures with 2 stems: 9555
    (overlaps: 0, !energy: 165582)
    </pre><P><H2>Sorting the structures...
    <P></H2><pre></pre><H2><P><P><P>Filtered and Sorted solutions:<P><P><P></H2><pre>
    </pre><H2><P><P><P><a HREF="http://biwww2.informatik.uni-freiburg.de/Software/MARNA/index.html" target="_blank">MARNA</a>-formatted:<P><P><P></H2><pre>
    GAACAUGGUUCUUGCCUUUUACCAGAACCAUCCGGGUGUUG
    ((((((((((..))))))))))((((((((...)))))))) -33.20 ( -0.69)
    (((((((((....)))))))))((((((((...)))))))) -33.17 ( -0.69)
    ((((((((((((((((...))))))))))))......)))) -32.40 ( +0.00)

"""
from __future__ import print_function

import sys
import os
import argparse
import pandas as pd
import progressbar
import logging


from rna_tools.tools.rna_alignment.rna_alignment import RNAalignment
from rna_tools.Seq import RNASequence

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def ungap(x):
    return x.replace('-', '')

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument('--one', help="one only for the first seq", action="store_true")
    parser.add_argument('--method', help="mcfold or rnastructure_CycleFold", default="mcfold")#rnastructure_CycleFold") #mcfold")
    parser.add_argument('--csv')
    parser.add_argument('--loop-seq', action="store_true")
    parser.add_argument('--template')
    parser.add_argument('--flanks', help="GC be default") #, default='GC')
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('alignment', help="an alignment in the Stockholm format")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    skipped = 0
    skipped_with_gaps = 0

    a = RNAalignment(args.alignment)

    df = pd.DataFrame()
    df_wide = pd.DataFrame()

    if args.template:
        f = open(args.template)
        lines = f.readlines()
        seqt = lines[0]
        sst = lines[1]
        cstl = lines[1:]  # a list of various cst

    if args.loop_seq:
        name = os.path.basename(args.alignment.replace('.stk', '').replace('.sto', '')) + \
          os.path.basename(args.template.replace('.cet', '')) + '-loop-seq'

    elif args.template:
        tp = os.path.basename(args.template.replace('.cet', ''))
        name = os.path.basename(args.alignment.replace('.stk', '').replace('.sto', ''))
        # + '-' + tp + '-upper-loop-'  + \               args.loop_upper + '-lower-loop-' + args.loop_lower
    if args.flanks:
        name += '-flanks-' + args.flanks

    try:
        name
    except:
        parser.print_usage()
        sys.exit(0)

    ################################################################################
    logging.basicConfig(level=logging.DEBUG,
                    format='%(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=name + '.log',
                    filemode='w')

    console = logging.StreamHandler()
    formatter = logging.Formatter('%(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


    log = logging.getLogger()
    if args.verbose:
        console.setLevel(logging.INFO)
    else:
        console.setLevel(logging.ERROR)
    log.setLevel(logging.INFO)

    log.info('name' + name)
    log.info(args)

    ################################################################################

    bar = progressbar.ProgressBar(max_value=len(a))
    bar.update(1)
    cc = 1
    for seq in a:
        cc += 1
        bar.update(cc)

        log.info('=' * 80)
        log.info(seq)
        log.info(seq.seq)
        # What to show? sst (from .cet file) of ss processed with RPT
        if args.template:
            log.info(sst)
        else:
            log.info(seq.ss)

        if 'N' in seq.seq: # + 'N':
            skipped += 1
            continue  # skip a seq with N

        method = args.method

        c = 1
        dfrow = ''
        for cst in cstl:
            seqq = RNASequence(seq.seq)
            # ech, remove gaps, but this RNASequence is not the same as RNASeq from the alignment
            seqq.seq = seqq.seq.replace('-', '')
            cst  = cst.replace('-', '').strip()
            if len(cst) != len(seqq.seq):
                continue
            energy, ss = seqq.predict_ss(method, constraints=cst, verbose=args.debug)
            log.info(ss + ' ' + str(energy) + ' ' + cst)
            energy_lower = energy

            df = df.append({
                'id': seq,
                'ss': ss,
                'cst': cst,
                'energy': energy,
                }, ignore_index=True)

            dfrow += ',"ss' + str(c) + '":"' + ss + '","energy' + str(c) + '":' + str(energy) + ',"cst' + str(c) + '":"' + cst + '"'
            c += 1

        dftxt = '{"id":"' + str(seq) + '"' + dfrow + '}'
        df_wide = df_wide.append(eval(dftxt), ignore_index=True)
        if args.one:
            if cc > 4:
                break

        # ok, anyway, save after each loop
        df.to_csv(name + '.csv')
        df_wide.to_csv(name + '_wide.csv')


    log.info(df)
    df.to_csv(name + '.csv')

    log.info(df_wide)
    df_wide.to_csv(name + '_wide.csv')

    log.info('Done: ' + name)
    print('Done: ' + name)
