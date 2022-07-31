#!/usr/bin/env python
"""mq scoring

A new column mqapRNA is created.

Input:

- <job_id>_raw.csv
- <job_id>_mqapRNA.csv
- <job_id>_raw_mqapRNA.csv

"""
from __future__ import print_function

import csv
import argparse
import pandas as pd
import os
import h2o

# plotting inside ipython
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

from sklearn.preprocessing import MinMaxScaler
from rna_tools.config import MODEL_KEY, MODEL_PATH
# plt.style.use('ggplot')#seaborn-deep')
import pandas as pd



def do_scoring(raw_csv):
    """Do scoring with h2o, reads output_csv and save to output to this csv.

    raw_csv of mq
    """

    print('in:', os.path.abspath(raw_csv))
    raw_csvf = os.path.abspath(raw_csv)
    model_key = MODEL_KEY
    print('model_key: ' + model_key)

    # why I need this?
    # path = os.getcwd()
    # if not path:
    #    path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(__file__)))
    # output_csv = path + os.sep + output_csv
    ##

    h2o.init(ip='localhost')

    # load the model
    h2o.load_model(
        MODEL_PATH
        )

    # load model into python
    df = h2o.import_file(path=raw_csvf)  # path + os.sep +
    print(df)
    dp_model = h2o.get_model(model_key)
    predict = dp_model.predict(df)
    df['mqapRNA'] = dp_model.predict(df)
    print(df['mqapRNA'])
    print(df.head()[0])

    # export it
    mq_csvf = raw_csvf.replace('.csv', '_mqapRNA.csv')  # <job_id>_raw_mqapRNA.csv
    print('Exporting to ', mq_csvf)
    h2o.export_file(df, mq_csvf, True)

    # sort it by mqapRNA
    df = pd.DataFrame.from_csv(mq_csvf)
    df_sort = df.sort_values(by=['mqapRNA'], ascending=[True])

    df_sort = df_sort.round(2)

    df[df['mqapRNA'] > 5] = 5
    # max score can be 5 (!!!!!!!!)
    ## scores = []
    ## for index, row in df.iterrows():
    ##     if row['mqapRNA'] > 5:
    ##          score = 5
    ##     else:
    ##          score = row['mqapRNA']
    ## # scores.append(score)
    ## df['mqapRNA'] = pd.Series(scores)

    print(df_sort)
    # save only mqapRNA column

    df_mini = pd.DataFrame(df_sort, columns=['mqapRNA'])
    df_mini.round(3)
    print(df_mini)
    out2 = raw_csvf.replace('_raw.csv', '_mqapRNA.csv')
    print('  saved ', out2)
    df_mini.to_csv(out2, sep="\t")

    h2o.cluster().shutdown()


def do_plot(raw_csvf, verbose=True):
    """

    Args:

       raw_csvf (str): /home/mqapRNA/mqapRNAweb/media/jobs/a49a2935-3571-44ea-9e0e-b0af5a469322/a49a293_raw.csv

    df = all scores
    """
    mq_csvf = raw_csvf.replace('.csv', '_mqapRNA.csv')  # <job_id>_raw_mqapRNA.csv
    df = pd.read_csv(mq_csvf)

    df.escore = -df.escore

    ## super ugly hack, if ss is error, then do ss_agree = 0
    try:
        df.ss_disagreement = -df.ss_disagreement
    except:
        df.ss_disagreement = 0
        pass
    ##

    df = df[['fn', 'ss_disagreement', 'clash_score', 'x3rnascore', 'simrna_total_energy', 'rnakb_lj_sr', 'escore',
             'farna_score_lowres', 'analyze_geometry', 'mqapRNA']]
    # normalize
    # cols that I want to use
    my_cols = set(df.columns)
    my_cols.remove('fn')
    my_cols = list(my_cols)
    df_nofn = df[my_cols]
    scaler = MinMaxScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df_nofn), columns=df_nofn.columns)

    df2 = pd.concat([df['fn'], df_scaled], axis=1)
    fns = []
    for index, row in df2.iterrows():
        fns.append(row['fn'][:20])
    df2['fn'] = pd.Series(fns)

    df2 = df2.sort_values(by=['mqapRNA'], ascending=[True])

    df2 =  df2.round(2)

    # try to plot thin lines and thing for mq
    #ax = df2[['fn', 'ss_disagreement', 'x3rnascore', 'simrna_total_energy', 'rnakb_lj_sr', 'escore',
    #          'farna_score_lowres', 'analyze_geometry']].plot(x='fn', ylim=(0, 1), rot=40, fontsize=10, linewidth=1)
    #
    ax = df2[['fn', 'mqapRNA']].plot(x='fn', ylim=(0, 1), rot=40, fontsize=10, linewidth=5)  # , ax=ax)
    ax.set_xlabel("Models")
    ax.set_ylabel("Normalized score. The lower the better quality.")

    plt.tight_layout()  # 'clash_score'
    plt.legend(loc=2, prop={'size': 8})

    plotfn = os.path.dirname(os.path.abspath(raw_csvf)) + os.sep + 'plot.png'
    if verbose:
        print(df2)
        print('ploting to %s [OK]' % plotfn)

    plt.savefig(plotfn, bbox_inches='tight', dpi=80)


# start
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='csv file (output mqaprna csv file)')
    args = parser.parse_args()
    do_scoring(args.input)
    # do_plot(args.input)
