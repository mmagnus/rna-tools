#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import pandas as pd
import argparse
import os

#plotting inside ipython
import matplotlib.pyplot as plt
import matplotlib
from sklearn.preprocessing import MinMaxScaler

matplotlib.style.use('ggplot')#seaborn-deep')
import pandas as pd

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('csv', help="", default="")
    # parser.add_argument('--sep', help="default is ,; can be also '\t'", default=",")
    # parser.add_argument("-v", "--verbose",
    #action="store_true", help="be verbose")
    parser.add_argument('--plot', help="output plot", default="plot.png")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    df = df.sort_values(by=['mqapRNArst'], ascending=[False])

    df['SS'] = df['mqapRNA'] + df['inf'] + df['rst_score']
    df['RST'] = df['mqapRNA'] + df['rst_score']

    ax = df.plot(x='fn', y="SS", kind="barh", color="orange")
    ax = df.plot(x='fn', y="RST", kind="barh", ax=ax, color="gray")#, figsize=(10,10));
    ax = df.plot(x='fn', y="mqapRNA", kind="barh", ax=ax, color="darkred")

    ax.set_title('mqapRNA score (the lower the better)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(args.plot, dpi=300)
