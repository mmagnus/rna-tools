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
from scipy.cluster import hierarchy
import numpy as np


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--font-scale', help="default=0.6", default=0.6)
    parser.add_argument('--sep', help="sep for input file, default=' '", default=",")
    parser.add_argument('--annot', help="annotate squares, default=False",
                        action="store_true", default=False)
    parser.add_argument('--plot', help="output plot, default=_plot_.png", default="_plot_.png")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    df = pd.read_csv(args.file,
                     sep=args.sep, index_col=False)#, index=False)
    print(df)
    mat = df.to_numpy()
    print(mat)
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
    Z = hierarchy.linkage(mat, 'single')
    matplotlib.rcParams['lines.linewidth'] = 2
    plt.style.use('dark_background')
    dn = hierarchy.dendrogram(Z)
    #ax.set_linecolor('orange')
    print(dn)
    plt.savefig(args.plot, dpi=300, transparent=True)
