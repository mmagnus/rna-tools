#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from rna_tools.tools.mq.RASP.RASP import RASP
from rna_tools.rna_tools_config import WRAPPERS_PATH, TEST_DATA


def main():
    wrapper = RASP('', '')
    try:
        result = wrapper.run(TEST_DATA + os.sep + '4ts2.pdb')
        if result:
            print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

if __name__ == '__main__':
    main()
