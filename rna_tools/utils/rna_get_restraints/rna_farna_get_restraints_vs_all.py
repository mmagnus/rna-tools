#!/usr/bin/python
import copy
import sys

import os
from subprocess import getoutput

f = sys.argv[1].strip().split(',')

"""
[ atompairs ]
C1' 12 C1' 61 FADE -25 25 10 -100 100
C1' 12 C1' 68 FADE -25 25 10 -100 100
C1' 61 C1' 68 FADE -25 25 10 -100 100
"""
print('[ atompairs ]')
f2 = copy.copy(f)
f2.reverse()

for i in f:
    for x in f2:
        if i !=  x:
            #print i,x
            print(("C1' %s C1' %s FADE -25 25 10 -100 100" % (i, x)))
    f2.pop()
