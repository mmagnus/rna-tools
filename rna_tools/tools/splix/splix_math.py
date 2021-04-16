#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse


"""
python splix_math.py 
[2, 1]
[0, 0]
[1, 0]
cwc15      1 1
prp16-302  1 0
bsc        -1 -1
bsg        0 -1
Δcwc15     1 1
"""

import gc

    
class Allele():
    def __init__(self, name, mechanism):
        self.name = name
        self.m = mechanism
        #if '\'' in mechanism:
        #    self.m = [1,0]
        #if '/' in mechanism:
        #    self.m = [0,1]

    def __add__(self, other):
        return [self.m[0] + other.m[0], self.m[1] + other.m[1]]  

    def __sub__(self, other):
        return [self.m[0] - other.m[0], self.m[1] - other.m[1]]  

    def delete(self):
        return Allele('Δ' + self.name, [self.m[1], self.m[0]])

    def __repr__(self):
        return ' '.join([self.name.ljust(10), str(self.m[0]), str(self.m[1])])
    #def __del__(self, other):
    #    return [-self.m[0], self.m[1]]


if __name__ == '__main__':
    # \ -> [0,5]
    c15 = Allele('cwc15',  [+1,+1])
    _302 = Allele('prp16-302',  [+1,0])
    bsc = Allele('bsc', [-1, -1])
    bsg = Allele('bsg', [0,-1])
    dc15 = c15.delete()

    print(dc15 + _302)
    print(c15 + bsc)
    print(c15 + bsg)

    for a in gc.get_objects(): # [c15, _302]:
        if isinstance(a, Allele):
            print(a)
