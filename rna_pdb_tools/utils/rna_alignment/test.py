#!/usr/bin/env python

"""quick and dirty test"""

import unittest
import subprocess
import os

class MyTest(unittest.TestCase):
    """RNA Alignment test"""
    def test_rna_align_find_seq_in_alignment(self):
        print('----------------------------------------------------------------------')
        cmd = "./rna_align_find_seq_in_alignment.py -a test_data/RF00167.stockholm.sto -f test_data/xx.fa"
        print(cmd)
        code = os.system(cmd)
        self.assertEqual(code, 0)
        
    def test_rna_align_find_core(self):
        print('----------------------------------------------------------------------')
        cmd = "./rna_align_find_core.py test_data/RF00167.stockholm.sto"
        print(cmd)
        code = os.system(cmd)
        self.assertEqual(code, 0)

    def test_rna_align_seq_to_alignment(self):
        print('----------------------------------------------------------------------')
        cmd = "./rna_align_seq_to_alignment.py -f test_data/4lvv_cmalign.txt -a test_data/RF01831.stockholm.sto"
        print(cmd)
        code = os.system(cmd)
        self.assertEqual(code, 0)

    def test_rna_align_seq_to_alignment_2(self):
        print('----------------------------------------------------------------------')
        cmd = "./rna_align_seq_to_alignment.py -s test_data/4lvv.seq -a test_data/RF01831.stockholm.sto -m test_data/RF01831.cm"
        print(cmd)
        code = os.system(cmd)
        self.assertEqual(code, 0)

    def test_rna_align_get_ss_from_alignment(self):
        print('----------------------------------------------------------------------')
        cmd = "./rna_align_get_ss_from_alignment.py test_data/ade.fa"
        print(cmd)
        code = os.system(cmd)
        self.assertEqual(code, 0)

    def test_rna_alignment(self):
        print('----------------------------------------------------------------------')
        cmd = "./rna_alignment.py"
        print(cmd)
        code = os.system(cmd)
        self.assertEqual(code, 0)
    def test_random_assignment_of_nucleotides(self):
        print('----------------------------------------------------------------------')
        cmd = "python random_assignment_of_nucleotides.py --alignfn test_data/aln1.fasta"
        print(cmd)
        code = os.system(cmd)
        self.assertEqual(code, 0)


if __name__ == '__main__':
    unittest.main()
