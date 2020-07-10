#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from PIL import Image, ImageChops, Image, ImageDraw, ImageFont, ImageStat
import logging
import os
import argparse
import matplotlib.pyplot as plt

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                            action="store_true", help="be verbose")
    parser.add_argument("-s", "--show",
                            action="store_true", help="be verbose")
    parser.add_argument("samples",
                        help="be verbose", type=int)
    parser.add_argument("conditions",
                        help="be verbose", type=int)

    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


def get_rms(im):
    stat = ImageStat.Stat(im)
    #r,g,b = stat.mean
    ## print('bg sum', stat.sum[0])
    ## print('bg mean', stat.mean[0])
    ## print('bg rms', stat.rms[0])
    return stat.rms[0]


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    img = Image.open(args.file)
    # design 8x6
    y = 0
    x = 0
    delx = img.width / args.samples
    dely = img.height / args.conditions

    dat = []
    for iy in range(0, args.conditions):
        vector = []
        for ix in range(0, args.samples):
            tmp = img.copy()
            if args.verbose: print((x, y, delx, dely))
            dot = tmp.crop((x, y, x+delx, y+dely))
            if args.show:
                dot.show()
            if ix == 0: # reset this is wt
                wt = get_rms(dot)
                vector.append(wt - wt + 0.5) # ugly hack to see anything
            else:
                sample = get_rms(dot)
                vector.append(wt - sample)           
            x += delx
        dat.append(vector)
        y += dely
        x = 0


    if args.verbose: print(dat)

    import numpy as np
    m = np.array(dat)
    if args.verbose: print(m)

    column_sums = m.sum(axis=0)
    print(column_sums)
    import pandas as pd
    df = pd.DataFrame({'score' : column_sums}, index = ['wt','del',34,39,40,410,41,42])
    df.plot.bar()

    from sklearn import preprocessing
    x = df.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled) #, index = 
    #df.set_index(['wt','del',34,39,40,410,41,42])
    p = print
    p(df)
    df.plot.bar()
    f = os.path.splitext(args.file)[0]
    outfn = f + '_plot.png'
    p('Output created ' + outfn)
    plt.savefig(outfn)
