#!/usr/bin/env python
r"""rna_plot_hist.py - generate a histogram

Don't open Excel, Jupyter. Simple plot a histogram of one column and save it to a file.

Example::

    # file
                                                      fn  rmsd_all
    0  19_Bujnicki_Human_4_rpr_n0-000001.pdb-000001_A...     14.73
    1  19_Bujnicki_Human_4_rpr_n0-000001.pdb.trafl_19...      0.46
    2  19_Bujnicki_Human_4_rpr_n0-000001.pdb.trafl_19...     14.73
    3  19_Bujnicki_Human_4_rpr_n0-000001.pdb_thrs0.50...      0.73
    4  19_Bujnicki_Human_4_rpr_n0-000001.pdb_thrs0.50...      0.83

    $ rna_plot_hist.py rmsds.csv --column rmsd_all

.. image:: ../../rna_tools/tools/plotting/test_data/rmsds_hist.png

"""
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import sys

plt.style.use('ggplot')
plt.rc('figure', figsize=(10, 6))
np.set_printoptions(precision=4, threshold=500)
pd.options.display.max_rows = 100


def get_parser():
    """Get parser."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="rmsd.txt")
    parser.add_argument('--column', help="column of file to plot")
    parser.add_argument('--sep', help="separator, be default \t", default=",")
    parser.add_argument('-o', '--output')
    parser.add_argument('--bins', default=10)
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    fn = args.file
    df = pd.read_csv(args.file, sep=args.sep)
    print(df.head())

    if not args.column:
        print('Set column')
        sys.exit(0)

    ax = df.plot.hist(args.column, grid=True, bins=int(args.bins))

    if not args.output:
        outfn = args.file.replace('.txt', '').replace('.csv', '') + '_hist.png'
        print ('Save plot %s' % outfn)
        plt.savefig(outfn)
