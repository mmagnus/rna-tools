#!/usr/bin/env python

"""
cmaling::

    [mm] thf cmalign RF01831.cm 4lvv.seq
    # STOCKHOLM 1.0
    #=GF AU Infernal 1.1.2

    4lvv         -GGAGAGUA-GAUGAUUCGCGUUAAGUGUGUGUGA-AUGGGAUGUCG-UCACACAACGAAGC---GAGA---GCGCGGUGAAUCAUU-GCAUCCGCUCCA
    #=GR 4lvv PP .********.******************9999998.***********.8999999******8...5555...8**************.************
    #=GC SS_cons (((((----(((((((((((,,,,,<<-<<<<<<<<___________>>>>>>>>>>,,,<<<<______>>>>,,,)))))))))))-------)))))
    #=GC RF      ggcaGAGUAGggugccgugcGUuAAGUGccggcgggAcGGGgaGUUGcccgccggACGAAgggcaaaauugcccGCGguacggcaccCGCAUcCgCugcc
    //

Reads seq files::
   
   >4lvv
   GGAGAGUAGAUGAUUCGCGUUAAGUGUGUGUGAAUGGGAUGUCGUCACACAACGAAGCGAGAGCGCGGUGAAUCAUUGCAUCCGCUCCA
"""
import rna_tools.utils.rna_alignment.rna_alignment as ra

import sys
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--file', help="cmalign output")
    parser.add_argument('-a', '--alignment', help="alignment file",  required=True)
    parser.add_argument('-m', '--cm', help="cm model to run cmalign on")
    parser.add_argument('-s', '--seq', help="seq fn, fasta!")
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args()

    # don't run cmalign
    if args.file:
        cma = ra.CMAlign(outputfn=args.file)
    else:
        cma = ra.CMAlign()#(outputfn=args.file)
        cma.run_cmalign(args.seq, args.cm)

    seq = cma.get_seq()

    a = ra.RNAalignment(args.alignment)
    print('cma hit  ' + seq)
    print('seq      ' + a.align_seq(seq))
    print('a.rf     ' + a.rf)



