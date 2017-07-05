#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Quick and ugly fix to update auto readme."""
from __future__ import print_function
import commands

if __name__ == '__main__':
    start_tag= "[mm] rna_pdb_tools$ git:(master) ✗ ./rna_pdb_tools.py -h"
    end_tag="  --delete DELETE"
    #
    doc = commands.getoutput('python rna_pdb_tools.py -h')
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
        if l.startswith(end_tag):
            copy_line_flag = True

    f = open('../README.md', 'w')
    f.write(fn)
    f.close()
