#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function


import argparse
import pathlib
import os
from PIL import Image, ImageChops, Image, ImageDraw, ImageFont, ImageStat, ImageOps


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--debug",
                        action="store_true")
    parser.add_argument("-a", "--levela", default = 50)
    parser.add_argument("-b", "--levelb", default = 60)
    parser.add_argument("-mx", default = 0, type=int)
    parser.add_argument("-my", default = 0, type=int)
    parser.add_argument("files", help="", default="", nargs='+')
    return parser

def trim(im):
    """
    https://stackoverflow.com/questions/10615901/trim-whitespace-using-pil
    """
    bg = Image.new(im.mode, im.size, im.getpixel((2,2)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 0.0, 0)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    else: return im


def intenisity():
    bg = Image.new(im.mode, im.size, im.getpixel((500,500)))
    diff = ImageChops.difference(im, bg)
    return diff


def level_image(img, a, b):
    f = '_tmp_.png'
    img.save(f)
    f2 = '_tmp_level_.png'
    cmd = 'convert ' + f + ' -level ' + str(a) + '%,' + str(b) +'% ' + f2 # -evaluate Min 90% out.png
    # print(cmd)
    os.system(cmd)
    img = Image.open(f2)

    os.remove(f)
    os.remove(f2)

    return img


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

    for f in args.files:
        print(f)

        rms = get_rms(Image.open(f))

        img = ImageOps.invert(Image.open(f))

        mask = Image.open(str(pathlib.Path(__file__).parent.absolute()) +  '/mask.png')
        #inverted_image.save('inv.png')

        img = img.convert('LA')
        if args.debug: img.save('_greyscale.png')

        img = img.transpose(Image.FLIP_LEFT_RIGHT)
        img = img.resize([1000, 1000])
        # img = normalized(img)
        img.paste(mask, (args.mx, args.my), mask)

        if rms < 180:
            img1 = level_image(img, args.levela, args.levelb)
            img1.save(os.path.splitext(f)[0] + '_prep.png')
        else:
            img2 = level_image(img, 40, 50)
            img2.save(os.path.splitext(f)[0] + '_prep.png')    
        #if args.debug: img.save('_before_trim.png')
        # img = trim(img)
        #img = img.resize([1000, 1000])
        #if args.debug: img.save('_after_trim.png')
        # if args.debug:
        #    img.save('_before_mask.png')
        # img.show()
