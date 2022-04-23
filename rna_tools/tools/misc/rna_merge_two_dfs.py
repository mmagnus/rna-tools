#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read more on ``--how`` (default: "inner")::

    how{‘left’, ‘right’, ‘outer’, ‘inner’}, default ‘inner’
    Type of merge to be performed.
    left: use only keys from left frame, similar to a SQL left outer join; preserve key order.
    right: use only keys from right frame, similar to a SQL right outer join; preserve key order.
    outer: use union of keys from both frames, similar to a SQL full outer join; sort keys lexicographically.
    inner: use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys.

Read more on optional ``--validate``::

    If specified, checks if merge is of specified type.
    “one_to_one” or “1:1”: check if merge keys are unique in both left and right datasets.
    “one_to_many” or “1:m”: check if merge keys are unique in left dataset.
    “many_to_one” or “m:1”: check if merge keys are unique in right dataset.
    “many_to_many” or “m:m”: allowed, but does not result in checks.

https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.merge.html

Example::

    $ rna_merge_dfs.py triple-db.csv triple-bk-contacts.csv triple
               triple  instances  clashes  near  exists
    0     cHH_cSH_AAA          0        0     0       0
    1     cHH_cSH_AAC          0        0     0       0
    2     cHH_cSH_AAG          0        0     0       0
    3     cHH_cSH_AAU          0        0     0       0
    4     cHH_cSH_ACA          0        0     0       0
    ...           ...        ...      ...   ...     ...
    6907  tWW_tsS_UGU          0        0     0       0
    6908  tWW_tsS_UUA          0        8     0       1
    6909  tWW_tsS_UUC          0        0     0       0
    6910  tWW_tsS_UUG          0        3     0       1
    6911  tWW_tsS_UUU          0        0     0       0

    [6912 rows x 5 columns]
                                              fn_bk  no_bk_contacts       triple
    0     Triple_cHH_cSH_AGA_rpr.pdb.NoBackbone.pdb               0  cHH_cSH_AGA
    1     Triple_cHH_cSH_AGC_rpr.pdb.NoBackbone.pdb               1  cHH_cSH_AGC
    2     Triple_cHH_cSH_AGG_rpr.pdb.NoBackbone.pdb               2  cHH_cSH_AGG
    3     Triple_cHH_cSH_AGU_rpr.pdb.NoBackbone.pdb               2  cHH_cSH_AGU
    4     Triple_cHH_cSH_CGA_rpr.pdb.NoBackbone.pdb               1  cHH_cSH_CGA
    ...                                         ...             ...          ...
    4188  Triple_tWW_tSW_UGU_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UGU
    4189  Triple_tWW_tSW_UUA_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UUA
    4190  Triple_tWW_tSW_UUC_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UUC
    4191  Triple_tWW_tSW_UUG_rpr.pdb.NoBackbone.pdb               4  tWW_tSW_UUG
    4192  Triple_tWW_tSW_UUU_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UUU

Example for outer::

    (py37) [mx] triples$ git:(master) ✗ rna_merge_dfs.py triple-db.csv triple-bk-contacts.csv triple --how outer
               triple  instances  clashes  near  exists
    0     cHH_cSH_AAA          0        0     0       0
    1     cHH_cSH_AAC          0        0     0       0
    2     cHH_cSH_AAG          0        0     0       0
    3     cHH_cSH_AAU          0        0     0       0
    4     cHH_cSH_ACA          0        0     0       0
    ...           ...        ...      ...   ...     ...
    6907  tWW_tsS_UGU          0        0     0       0
    6908  tWW_tsS_UUA          0        8     0       1
    6909  tWW_tsS_UUC          0        0     0       0
    6910  tWW_tsS_UUG          0        3     0       1
    6911  tWW_tsS_UUU          0        0     0       0

    [6912 rows x 5 columns]
                                              fn_bk  no_bk_contacts       triple
    0     Triple_cHH_cSH_AGA_rpr.pdb.NoBackbone.pdb               0  cHH_cSH_AGA
    1     Triple_cHH_cSH_AGC_rpr.pdb.NoBackbone.pdb               1  cHH_cSH_AGC
    2     Triple_cHH_cSH_AGG_rpr.pdb.NoBackbone.pdb               2  cHH_cSH_AGG
    3     Triple_cHH_cSH_AGU_rpr.pdb.NoBackbone.pdb               2  cHH_cSH_AGU
    4     Triple_cHH_cSH_CGA_rpr.pdb.NoBackbone.pdb               1  cHH_cSH_CGA
    ...                                         ...             ...          ...
    4188  Triple_tWW_tSW_UGU_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UGU
    4189  Triple_tWW_tSW_UUA_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UUA
    4190  Triple_tWW_tSW_UUC_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UUC
    4191  Triple_tWW_tSW_UUG_rpr.pdb.NoBackbone.pdb               4  tWW_tSW_UUG
    4192  Triple_tWW_tSW_UUU_rpr.pdb.NoBackbone.pdb               3  tWW_tSW_UUU

    [4193 rows x 3 columns]
                        triple  instances  clashes  near  exists                                              fn_bk  no_bk_contacts
    0              cHH_cSH_AAA        0.0      0.0   0.0     0.0                                                NaN             NaN
    1              cHH_cSH_AAC        0.0      0.0   0.0     0.0                                                NaN             NaN
    2              cHH_cSH_AAG        0.0      0.0   0.0     0.0                                                NaN             NaN
    3              cHH_cSH_AAU        0.0      0.0   0.0     0.0                                                NaN             NaN
    4              cHH_cSH_ACA        0.0      0.0   0.0     0.0                                                NaN             NaN
    ...                    ...        ...      ...   ...     ...                                                ...             ...
    7228           tWW_tSS_UGA        NaN      NaN   NaN     NaN          Triple_tWW_tsS_UGA_rpr.pdb.NoBackbone.pdb             3.0
    7229           tWW_tSS_UGG        NaN      NaN   NaN     NaN          Triple_tWW_tsS_UGG_rpr.pdb.NoBackbone.pdb             4.0
    7230           tWW_tSS_UUA        NaN      NaN   NaN     NaN          Triple_tWW_tsS_UUA_rpr.pdb.NoBackbone.pdb             0.0
    7231           tWW_tSS_UUG        NaN      NaN   NaN     NaN          Triple_tWW_tsS_UUG_rpr.pdb.NoBackbone.pdb             0.0
    7232  tWW_tSW_CGC_exemplar        NaN      NaN   NaN     NaN  Triple_tWW_tSW_CGC_exemplar_rpr.pdb.NoBackbone...             7.0

"""
from __future__ import print_function
import pandas as pd
import argparse
import os

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('csv1', help="", default="")
    parser.add_argument('--sep1', help="default is ,; can be also '\t'", default=",")
    parser.add_argument('csv2', help="", default="")
    parser.add_argument('--sep2', help="default is ,; can be also '\t'", default=",")
    parser.add_argument('mergeon', help="merge on column", default="")
    parser.add_argument('--how', help="inner | left | right | outer, default inner", default="inner")
    parser.add_argument('--validate', help="default: one_to_one", default="one_to_one")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument('--output', help="output csv", default="merged.csv")
    parser.add_argument('--force-writing-output',  action="store_true")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    df1 = pd.read_csv(args.csv1, delimiter=args.sep1)
    try:
        df1.pop('id')
    except:
        pass
    df2 = pd.read_csv(args.csv2, delimiter=args.sep2)
    try:
        df2.pop('id')
    except:
        pass
    print(df1)
    print(df2)

    if args.validate:
        merged = pd.merge(df1, df2, on=args.mergeon, how=args.how, validate=args.validate)
    else:
        merged = pd.merge(df1, df2, on=args.mergeon, how=args.how)
    print(merged)
    if not os.path.isfile(args.output) or args.force_writing_output:
        merged.to_csv(args.output, index=False)
        print('Saved to %s' % args.output)
    else:
        print('! Output file exist, NOT SAVED (remove the file and run again) %s' % args.output)
