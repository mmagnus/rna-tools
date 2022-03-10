#!/usr/bin/python
# -*- coding: utf-8 -*-
from unittest import TestCase, main
from clarna_wrapper import ClaRNAWrapper
import os

class Tests(TestCase):
    def setUp(self):
        os.system('./doc_generate_ex_output.sh')

    def test_Basic_1xjr(self):
        f = 'test_data/1mme.pdb'

        c = ClaRNAWrapper(f)
        seq = c.get_seq()
        ss = c.get_ss()

        self.assertEqual(seq, 'UGGUGCUGAUGAGGCCACCACGUAACCGGGAAACUCUGAGCUGGUGCUGAUGAGGCCACCACGUAACCGGGAAACUCUGAGC')
        self.assertEqual(ss, '(((.(.......(((())).)....))))....(((.))).(((((.......(((()))))....))))....(((.))).')

if __name__ == '__main__':
    main()
