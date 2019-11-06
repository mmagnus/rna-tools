#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
rna_calc_rmsd_multi_targets.py - calculate RMSDs of structures to multiple targets::

   $ rna_calc_rmsd_multi_targets.py --models multi-targets/rp21/*.pdb 
                                    --targets multi-targets/rp21/solutions/*.pdb 
                                    --target-selection A:1-27+29-41 
                                    --model-selection A:1-27+29-41

CSV table produced::

                              21_solution_0_ChainA.pdb  21_solution_0_ChainB.pdb  21_solution_1_ChainA.pdb  21_solution_1_ChainB.pdb  21_solution_2.pdb   mean    min    max    sd
    fn
    21_3dRNA_1_rpr.pdb                           12.17                     12.11                     12.17                     12.11              12.11  12.13  12.11  12.17  0.03
    21_Adamiak_1_rpr.pdb                          4.64                      4.61                      4.64                      4.61               4.64   4.63   4.61   4.64  0.01
    21_ChenHighLig_1_rpr.pdb                      4.01                      3.97                      4.01                      3.97               4.07   4.01   3.97   4.07  0.04
    21_Das_1_rpr.pdb                              5.71                      5.60                      5.71                      5.60               5.61   5.65   5.60   5.71  0.05
    --------------------------------------------------------------------------------
    Save rna_calc_rmsd_multi_targets_output.csv


"""
import argparse
import os
import pandas as pd


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # parser.add_argument('-args', help="", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--models", help="", nargs='+')
    parser.add_argument("--targets", help="", nargs='+')
    parser.add_argument("--output-csv", help="", default='rna_calc_rmsd_multi_targets_output.csv')
    parser.add_argument("--model-selection",
                         default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")
    parser.add_argument("--target-selection",
                         default='',
                         help="selection, e.g. A:10-16+20, where #16 residue is included")
    
    return parser


def hr(t):
    print
    if t:
        print(t)
    print('-' * 80)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    
    if args.verbose:
        print(args.models)
        print(args.targets)

    rmsds_files = []
    # run rmsd_calc_rmsd.py
    for target in args.targets:
        hr(target)
        outfile = os.path.basename(target).replace('.pdb', '.csv')
        cmd = "rna_calc_rmsd.py -t " + target + " " + ' '.join(args.models) + ' -o ' + outfile + ' --target-column-name '
        if args.target_selection:
            cmd += ' --target-selection ' + args.target_selection
        if args.model_selection:
            cmd += ' --model-selection ' + args.model_selection
        if args.verbose: print(cmd)
        os.system(cmd)
        rmsds_files.append(outfile)
    # process the results
    dfs = []
    for f in rmsds_files:
        df1 = pd.read_csv(f, index_col='fn')
        dfs.append(df1)

    df = pd.concat(dfs, axis=1)# join='inner')
    mean = df.mean(numeric_only=True, axis=1)
    median = df.median(numeric_only=True, axis=1)
    var = df.var(numeric_only=True, axis=1)
    min = df.min(numeric_only=True, axis=1)
    max = df.max(numeric_only=True, axis=1)
    sd = df.std(numeric_only=True, axis=1)

    df['mean'] = mean
    df['median'] = median
    df['variance'] = var
    df['min'] = min
    df['max'] = max
    df['sd'] = sd
    
    df.round(2).to_csv(args.output_csv)
    print('-' * 80)
    print(df.round(2))
    print('-' * 80)
    print('Save %s' % args.output_csv)
