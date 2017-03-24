#!/usr/bin/python
#-*- coding: utf-8 -*-

import sys
import numpy
import optparse

chars = ['.', '▁', '▂', '▃', '▄', '▅', '▆', '▇', '█']

if __name__ == '__main__':
    optparser=optparse.OptionParser(usage="%prog [<options>]", description="dupa", version="beta")
    (opts, args)=optparser.parse_args()

    if len(sys.argv) != 2:
        print((optparser.format_help()))
        sys.exit(1)


    fn = sys.argv[1]
    f = open(fn)

    # collect reactivites
    reactivites = []
    for l in f:
        reactivites.extend([float(x) for x in l.split()])

    freq, bins = numpy.histogram(reactivites, 7)
    bins = numpy.digitize(reactivites, bins)

    plot = ''
    for b in bins:
        # print chars[b],
        plot += chars[b]

    print(plot)