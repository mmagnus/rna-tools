#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Quick and ugly fix to update auto readme."""

from __future__ import print_function
from future import standard_library
import subprocess
standard_library.install_aliases()

if __name__ == '__main__':
    start_tag = "[mm] rna_pdb_tools$ git:(master) ✗ ./rna_pdb_toolsx.py -h"
    end_tag = "more than one fragment --extract"
    #
    doc = subprocess.getoutput('python rna_pdb_toolsx.py -h')
    print(doc)
    fn = ''
    f = open('../README.md', 'r')
    copy_line_flag = True
    for l in f:
        if copy_line_flag:
            fn += l
        if l.startswith(start_tag):
            fn += doc + '\n'
            copy_line_flag = False
        if l.strip().startswith(end_tag):
            copy_line_flag = True

    f = open('../README.md', 'w')
    f.write(fn)
    f.close()
