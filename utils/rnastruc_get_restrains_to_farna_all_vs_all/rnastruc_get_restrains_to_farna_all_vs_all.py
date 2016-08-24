#!/usr/bin/python
import copy
import sys

import os
from commands import getoutput
def get_version(currfn=__file__, verbose=False):
    """Get version of the tool based on state of the git repository.
    Return version. 
    If currfn is empty, then the path is '.'. Hmm.. I think it will work. We will see.
    The version is not printed!"""
    if currfn == '':
        path = '.'
    else:
        path = os.path.dirname(currfn)
    if verbose: print 'get_version::path', path
    if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
        path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))
    if not path: path = '.'
    if verbose: print 'get_version::path2', path
    curr_path = os.getcwd()
    os.chdir(os.path.abspath(path))
    version =getoutput('git describe --long --tags --dirty --always')
    os.chdir(curr_path)
    return version

f = sys.argv[1].strip().split(',')

"""
[ atompairs ]
C1' 12 C1' 61 FADE -25 25 10 -100 100
C1' 12 C1' 68 FADE -25 25 10 -100 100
C1' 61 C1' 68 FADE -25 25 10 -100 100
"""
print get_version()
print '[ atompairs ]'
f2 = copy.copy(f)
f2.reverse()

for i in f:
    for x in f2:
        if i !=  x:
            #print i,x
            print "C1' %s C1' %s FADE -25 25 10 -100 100" % (i, x)
    f2.pop()
