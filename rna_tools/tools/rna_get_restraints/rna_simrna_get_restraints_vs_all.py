#!/usr/bin/python
import copy
import sys

import os
from subprocess import getoutput

f = sys.argv[1].strip().split(',')

f2 = copy.copy(f)
f2.reverse()

txt = ''
for i in f:
    for x in f2:
        if i !=  x:
            #print i,x
            #print "C1' %s C1' %s FADE -25 25 10 -100 100" % (i, x)
            txt += 'SLOPE A/' + str(i) + '/MB A/' + str(x) + '/MB 15 25 0.5 \n' # energy values about restrains. slope <-> penalty
            txt += 'WELL A/' + str(i) + '/MB A/' + str(x) + '/MB  15 25 0.5 \n' # well <-> reward # check paper
    f2.pop()
print((txt.strip()))
