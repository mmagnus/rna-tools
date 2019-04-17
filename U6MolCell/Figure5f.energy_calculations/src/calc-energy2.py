#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate energy (.cet) format::

    UGGC-CCCUGCGCAA-GGAUGACA
    (((..((((.....).))..))))
    lll---uuuoooooU-UU---LLL
    (((..(((((***)).))..))))

Examples::

    calc-energy2.py --template alignments/u6-lower.cet alignments/u6-only-RemovedGapped.stk -v
          --loop-upper guaa --loop-lower guaa
          --loop-upper-cst '(..)' --loop-lower-cst '(..)'
    calc-energy2.py --template u6atac-template.txt u6atac_u6only.sto -v
   ./calc-energy2.py --template alignments/u6-lower.cet --one alignments/u6-lower-stem-only.sto

Takes cet files (calc-energy-templets).

    [mm] energy calc box$ git:(master) ✗ ./calc-energy2.py --template alignments/u6-lower.cet --one alignments/u6-only.stk  -v --loop-seq
    u6-only-loop-seq-u6-lower
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


"""
from __future__ import print_function
import pandas as pd
import sys
import argparse
import progressbar
import os
import logging

import pandas as pd
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
    #parser.add_argument('--csv', help="", default="tmp.csv")
    parser.add_argument('--method', help="mcfold or rnastructure_CycleFold", default="mcfold")#rnastructure_CycleFold") #mcfold")
    parser.add_argument('--model', help="4x4, 3x3")#, action="store_true")
    #parser.add_argument('--m3x3', help="3x3 model", action="store_true")
    parser.add_argument('--plot-each-step', action="store_true")
    parser.add_argument('--loop-upper', default="guaa") # cau")
    parser.add_argument('--loop-lower', default="guaa")
    parser.add_argument('--loop-upper-cst', default="(**)")
    parser.add_argument('--loop-lower-cst', default="(**)")
    parser.add_argument('--csv')
    parser.add_argument('--loop-seq', action="store_true")
    parser.add_argument('--template')
    parser.add_argument('--flanks', help="GC be default") #, default='GC')
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('--stop-if-size-diff', action="store_true")
    #parser.add_argument('--name', help="file plot name (+.png), csv name (+.csv), title of the plot", default="tmp")
    #parser.add_argument('--title')
    #, default='U6: energies of lower vs upper stem')
    parser.add_argument('alignment', help="an alignment in the Stockholm format")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    from rna_tools.tools.rna_alignment.rna_alignment import RNAalignment
    from rna_tools.Seq import RNASequence

    skipped = 0
    skipped_with_gaps = 0

    a = RNAalignment(args.alignment)

    df = pd.DataFrame()

    if args.template:
        f = open(args.template)
        seqt = f.readline()  # t - from template
        sst = f.readline().strip()  # sst from template
        template = f.readline()
        cstl = f.readline()
        # ---tttllll-uuuu-ooooooooo-uuuu---llllttt--------
        #template = '---tttllll-uuuu-ooooooooo-uuuu---LLLLttt--------'

    ################################################################################
    if args.loop_seq:
        name = os.path.basename(args.alignment.replace('.stk', '').replace('.sto', '')) + \
          os.path.basename(args.template.replace('.cet', '')) + '-loop-seq'

    elif args.template:
        tp = os.path.basename(args.template.replace('.cet', ''))
        name = os.path.basename(args.alignment.replace('.stk', '').replace('.sto', '')) + '-' + tp + '-upper-loop-'  + \
               args.loop_upper + '-lower-loop-' + args.loop_lower
    if args.flanks:
        name += '-flanks-' + args.flanks
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

    log.info(name + '<= name')
    log.info(args)

    ################################################################################

    bar = progressbar.ProgressBar(max_value=len(a))
    bar.update(0)
    c = 0

    for seq in a:
        c += 1
        bar.update(c)

        log.info('=' * 80)
        log.info(seq)
        log.info(seq.seq)
        # What to show? sst (from .cet file) of ss processed with RPT
        if args.template:
            log.info(sst)
            log.info(template)
        else:
            log.info(seq.ss)

        if 'N' in seq.seq: # + 'N':
            skipped += 1
            continue  # skip a seq with N

        # loop seq
        # take a loop from the seq and use it for both, lower and upper stems
        loop_seq = seq.seq[9:14].lower()

        method = args.method
        # LOWER #
        log.info('lower ------------------------------------')
        if args.loop_seq:
            loop_lower = loop_seq # from alignment
        else:
            loop_lower = args.loop_lower

        if args.template:
            a = ''
            b = ''
            loop_lower = ''
            for s, t in zip(seq.seq, template):
                if t == 'l' and s != '-' :
                    a += s
                if t == 'L' and s != '-':
                    b += s
                if t == 'o' and s != '-':
                    loop_lower += s.lower()

            if not args.loop_seq:
                loop_lower = args.loop_lower

            if not args.loop_seq:
                loop_lower = args.loop_lower

            cst_a = ''
            cst_b = ''
            cst_l = ''
            for s, t in zip(cstl, template):
                if t == 'l' and s != '-' :
                    cst_a += s
                if t == 'L' and s != '-':
                    cst_b += s
                if t == 'o' and s != '-':
                    cst_l += s.lower()

            if args.loop_lower_cst:
                cst_l = args.loop_lower_cst

            cst = cst_a + cst_l + cst_b

        if args.flanks:
            a = args.flanks[0] + a
            b += args.flanks[1]
            cst = '(' + cst + ')'

        log.info(a + ' ' + loop_lower + ' ' + b)
        log.info('')
        log.info(a + loop_lower + b)
        log.info(cst + ' <= cst')

        # log.info('\_ loop lower: %s' % loop_lower)

        seql = RNASequence(a + loop_lower + b)
        #if '-' in seqn.seq: # + 'N':
        #    skipped_with_gaps += 1
        #    continue  # skip a seq with N
        #log.info(seql)
        energy, ss = seql.predict_ss(method, constraints=cst, verbose=args.debug)
        log.info( seql.seq + '\n' + ss + ' ' + str(energy))
        energy_lower = energy

        log.info('upper ------------------------------------')
        ## if args.loop_seq:
        ##     loop_upper = seq.seq[9:14].lower() # # from alignment
        ##     assert loop_upper, 'Loop is empty'
        ## else:
        ##     loop_upper = args.loop_upper
        ## log.info('dddd')
        ## log.info(loop_upper)
        if args.template:
            a = ''
            b = ''
            loop_upper = ''
            log.info(seq.seq)
            for s, t in zip(seq.seq, template):
                if t == 'u' and s != '-' :
                    a += s
                if t == 'U' and s != '-':
                    b += s
                if t == 'o' and s != '-':
                    loop_upper += s.lower()

            if not args.loop_seq:
                loop_upper = args.loop_upper

            cst_a = ''
            cst_b = ''
            cst_l = ''
            for s, t in zip(cstl, template):
                if t == 'u' and s != '-' :
                    cst_a += s
                if t == 'U' and s != '-':
                    cst_b += s
                if t == 'o' and s != '-':
                    cst_l += s.lower()

            if args.loop_upper_cst:
                cst_l = args.loop_upper_cst

            cst = cst_a + cst_l + cst_b

        if args.flanks:
            a = args.flanks[0] + a
            b += args.flanks[1]
            cst = '(' + cst + ')'
        # log.info('\_ loop lower: %s' % loop_lower)

        assert loop_upper, 'Loop is empty'
        log.info(a + ' ' + loop_upper + ' ' + b)
        log.info('')
        log.info(a + loop_upper + b)
        log.info(cst + ' <= cst')

        # upper
        seqn = RNASequence(a + loop_upper + b)

        # log.info('\_ loop upper: %s' % loop_upper)
        #log.info(seqn)
        energy_upper, ss_upper = seqn.predict_ss(method, constraints=cst, verbose=args.debug)
        log.info(seqn.seq + '\n' + ss_upper + ' ' + str(energy_upper))

        df = df.append({
            'id': seq.id,
            'low_seq': seql.seq,
            'low_ss': ss,
            'low_energy': energy_lower,
            'up_seq': seqn.seq,
            'up_ss': ss_upper,
            'up_energy': energy_upper,
            }, ignore_index=True)

        if args.plot_each_step:
            ax = df.boxplot()
            ax.set_title(name)
            fig = ax.get_figure()
            fig.savefig(name + '.png')
            df.to_csv(name + '.csv')

        if args.one:
            break

        if len(seql.seq) != len(seqn.seq):
            # raise Exception('Seqs of different size')
            log.info('Seqs of different size')
            if args.stop_if_size_diff:
                sys.exit(1)
            continue

    log.info(df)
    ax = df.boxplot()
    ax.set_title(name)
    fig = ax.get_figure()
    fig.savefig(name + '.png')
    df.to_csv(name + '.csv')
    log.info('Done: ' + name)
    print('Done: ' + name)
    #log.info('loop:' + loop + 'flank: ' + flank + 'common flanks:' + cflank)
    #log.info('Skipped', skipped)
    #log.info('skipped_with_gaps', skipped_with_gaps)
