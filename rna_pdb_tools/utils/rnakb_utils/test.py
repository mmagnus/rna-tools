#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .rnakb_utils import *

def test_gromacs_ready(fn):
    """Test if gromacs ready
    """
    cmd = 'pdb2gmx -ff amber03 -water none -f %s -o query.gro -p query.top' % fn
    os.system(cmd)

def test():
    """Test #1
    """
    fn = 'test_data/1msy_output4_01-000001_AA.pdb'
    f = open(fn).read()
    out = make_rna_rnakb_ready(f)
    o = open('test_output/1msy_output4_01-000001_AA_rnakb_ready.pdb', 'w')
    o.write(out)
    o.close()
    #test_gromacs_ready('tmp.pdb')

def test2():
    """Test #2
    """
    fn = 'test_data/1duq.pdb'
    f = open(fn).read()
    out = make_rna_rnakb_ready(f)
    o = open('test_output/1duq_rnakb_ready.pdb', 'w')
    o.write(out)
    o.close()

def test2b():
    fn = 'test_data/1duq_rpr.pdb'
    f = open(fn).read()
    out = make_rna_rnakb_ready(f)
    o = open('test_output/1duq_rpr_rnakb_ready.pdb', 'w')
    o.write(out)
    o.close()
    #test_gromacs_ready(fn)

def test3():
    """Test #3
    """
    fn = 'test_data/cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_error_fx.pdb'
    f = open(fn).read()
    out = make_rna_rnakb_ready(f)
    o = open('test_output/cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_error_fx_rnakb_ready.pdb', 'w')
    o.write(out)
    o.close()
    #test_gromacs_ready('tmp3.pdb')

def test4():
    """Test #3
    """
    fn = 'test_data/3nt_edited.pdb'
    f = open(fn).read()
    out = make_rna_rnakb_ready(f)
    o = open('test_output/3nt_edited_rnakb_ready.pdb', 'w')
    o.write(out)
    o.close()
    #test_gromacs_ready('tmp3.pdb')

if __name__ == '__main__':
    test()
    test2()
    test2b()
    test3()
    test4()
