#!/usr/bin/python

from pdb_checker import *
import sys

if __name__ == '__main__':
    fn = sys.argv[1]
    print 'Input: ', fn
    #fn = 'test_data/decoy0165_amb.pdb'
    t  = sys.argv[1].replace('.pdb', '_clx.pdb') 

    fix_rresnumes(fn, t)
    remove_hydrogen(t,t)
    remove_ion(t,t)
    remove_water(t,t)
    renum_atoms(t,t)
    fix_O_in_UC(t,t)
    fix_op_atoms(t, t)
