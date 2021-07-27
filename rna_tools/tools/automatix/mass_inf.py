#!/usr/bin/env python
import os
import glob

def main():
    root = os.getcwd()
    pdbs = glob.glob('*')
    for p in pdbs:
        if os.path.isdir(p):
            os.chdir(p) # go inside a folder
            cmd = 'rna_calc_inf.py -t struc/' + p + '_M1.pdb struc/*.pdb'
            print cmd
            #os.system
            os.chdir(root)

if __name__ == '__main__':
    main()
