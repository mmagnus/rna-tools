#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Quick and ugly fix to update auto readme."""

import subprocess

if __name__ == '__main__':
    start_tag = "usage: rna_pdb_toolsx.py [-h] [--version] [-r] [--renum-atoms]"
    end_tag = "ACCCGCAAGGCCGACGGC GCCGCCGCUGGUGCAAGUCCAGCCACGCUUCGGCGUGGGCGCUCAUGGGU"
    #
    doc = subprocess.Popen('python rna_pdb_toolsx.py -h', shell=True, stdout=subprocess.PIPE).stdout.read().decode()
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
