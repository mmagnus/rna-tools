#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from wrappers.ClashScore.ClashScore import ClashScore

def main():
    wrapper = ClashScore('', '')
    fns = ['test' + os.sep + '1xjrA.pdb', 
           'test' + os.sep + '1xjrA_M1.pdb',  # native, and M1 from rasp decoys
           'test' + os.sep + '6TNA.pdb'
    ]
    for f in fns:
        result = wrapper.run(f)
        print(f)
        print(result)
        wrapper.cleanup()

if __name__ == '__main__':
    main()
