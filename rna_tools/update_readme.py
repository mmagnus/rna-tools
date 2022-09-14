#!/usr/bin/env python3
"""Quick and ugly fix to update auto readme."""

import subprocess

if __name__ == '__main__':
    start_tag = "usage: rna_pdb_tools.py"
    end_tag = "Tricks:"

    #
    doc = subprocess.Popen('python rna_pdb_tools.py -h', shell=True, stdout=subprocess.PIPE).stdout.read().decode()
    txt = ''
    f = open('../README.md', 'r')
    copy_line_flag = True
    for l in f:
        if l.strip().startswith(end_tag):
            copy_line_flag = True
        if l.startswith(start_tag):
            txt += doc + '```\n\n'
            copy_line_flag = False
        if copy_line_flag:
            txt += l

    f = open('../README.md', 'w')
    f.write(txt)
    print(txt)
    f.close()
