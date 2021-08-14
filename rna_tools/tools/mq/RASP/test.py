#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from rna_tools.tools.mq.RASP.RASP import RASP

def main():
    wrapper = RASP()
    try:
        result = wrapper.run('../../../input/4ts2.pdb')
        if result:
            print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

if __name__ == '__main__':
    main()
