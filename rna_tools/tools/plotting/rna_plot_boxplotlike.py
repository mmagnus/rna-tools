#!/usr/bin/env python
r"""rna_plot_density.py - generate a density plot

Don't open Excel, Jupyter. Simple plot a density of one column and save it to a file.

Example::

    # file
                                                      fn  rmsd_all
    0  19_Bujnicki_Human_4_rpr_n0-000001.pdb-000001_A...     14.73
    1  19_Bujnicki_Human_4_rpr_n0-000001.pdb.trafl_19...      0.46
    2  19_Bujnicki_Human_4_rpr_n0-000001.pdb.trafl_19...     14.73
    3  19_Bujnicki_Human_4_rpr_n0-000001.pdb_thrs0.50...      0.73
    4  19_Bujnicki_Human_4_rpr_n0-000001.pdb_thrs0.50...      0.83

    $ rna_plot_hist.py rmsds.csv --column rmsd_all

.. image:: ../../rna_tools/tools/plotting/test_data/rmsds_dens.png
"""
from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys

plt.style.use('ggplot')
plt.rc('figure', figsize=(10, 6))


def get_parser():
    """Get parser."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help="rmsd.txt")
    parser.add_argument('x', help="column of file to plot")
    parser.add_argument('y', help="column of file to plot")
    parser.add_argument('--sep', help="separator, be default \t", default=",")
    parser.add_argument('-o', '--output')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    fn = args.file
    x = args.x
    y = args.y
    df = pd.read_csv(args.file, sep=args.sep)
    print(df.head())

    import pandas as pd
    import seaborn as sns 
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6,4))
    sns.set(style="white")
    #ax = sns.violinplot(data=df, x=x, y=y, color="orange")
    ax = sns.boxplot(data=df, x="growthb", y="rmsd_all", color="black")#, ax=ax)
    #ax.set_ylim(0,3)
    ax = sns.swarmplot(data=df, x="growthb", y="rmsd_all", color="orange", ax=ax)
    #plt.savefig('/Users/magnus/d/out.png', dpi=200)
    plt.tight_layout()
    if not args.output:
        outfn = args.file.replace('.txt', '').replace('.csv', '') + '_bx.png'
        print('Save plot %s' % outfn)
        plt.savefig(outfn, dpi=100)
    import os
    os.system('open %s' % outfn)
