#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""rna_cal_rmsd_trafl_plot - generate a plot based of <rmsd.txt> of rna_calc_evo_rmsd.py."""

from __future__ import division

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pandas import Series, DataFrame

plt.style.use('ggplot')
plt.rc('figure', figsize=(10, 6))
np.set_printoptions(precision=4, threshold=500)
pd.options.display.max_rows = 100


def get_parser():
    """Get parser."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="rmsd.txt")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    fn = args.file
    df = pd.read_csv(fn, sep="\t")
    print(df.head())
    print(len(df))  # print len(df) #

    # plt.style.use('classic')
    ax = df.plot(x=df.columns[0], y=df.columns[1],
                 legend=False,
                 marker="o",
                 clip_on=False)  # (x='fn', y='rmsd')
    ax.set_ylabel(df.columns[1])
    plt.ylim(-0, df[df.columns[1]].max() + 1)
    plt.xlim(-0, df[df.columns[0]].max() + 1)

    outfn = args.file.replace('.txt', '') + '.png'
    print('Save plot %s' % outfn)
    plt.savefig(outfn)
