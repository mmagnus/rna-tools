#!/usr/bin/python
# -*- coding: utf-8 -*-
from unittest import TestCase, main
from rna_x3dna import x3DNA, x3DNAMissingFile

class Tests(TestCase):
    def setUp(self):
        pass

    def test_Basic_1xjr(self):
        p = x3DNA('test_data/1xjr.pdb')
        seq = p.get_seq()
        self.assertEqual(seq,
                         'gGAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU'
                         )

    def test_FilesDoesNotExists(self):
        self.assertRaises(x3DNAMissingFile, lambda : \
                          x3DNA('test/dupa666'))


if __name__ == '__main__':
    main()
