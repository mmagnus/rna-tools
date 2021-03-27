#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
import matplotlib.pyplot as plt  # to start
import matplotlib

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    matplotlib.rcParams['lines.linewidth'] = 10
    plt.figure(figsize=(3,3))
    df = pd.read_csv(args.file,
                     sep=' ', index_col=False)#, index=False)
    df = df.set_index(df.columns)
    print(df)
    #df = df.drop('Unnamed: 15', axis=1)
    g = sns.clustermap(df, cmap="bwr", annot=True)#, linecolor = 'black', 
                           #linewidths=.3)#, vmin=-6, vmax=+6) # , 
    for a in g.ax_row_dendrogram.collections:
        a.set_linewidth(1)
        a.set_color('orange')
    for a in g.ax_col_dendrogram.collections:
        a.set_linewidth(1)
        a.set_color('orange')
    plt.savefig('_tmp.png', dpi=300)
