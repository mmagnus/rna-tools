#!/usr/bin/env python
#-*-coding: utf-8-*-

"""
tests for utils.py.
"""

__author__ = "Tomasz Puton"
__credits__ = "Kristian Rother"
__license__ = "GNU GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

import tempfile, os
from unittest import TestCase,  main

from .secstruc import BasePairs
from .parsers import parse_bp
from .writers import BpseqWriter, BpseqWriterError

class BpseqWriterTests(TestCase):
    def setUp(self):
        """Sets up environment for tests"""
        self.fd, self.path = tempfile.mkstemp()
        self.bplist1 = BasePairs([(2,5),(1,4),(9,13)])
        self.seq1 = 'ACACACACACACA'
        self.ref_bpseq1 = \
"""1 A 4
2 C 5
3 A 0
4 C 1
5 A 2
6 C 0
7 A 0
8 C 0
9 A 13
10 C 0
11 A 0
12 C 0
13 A 9
"""

        self.bplist2 = BasePairs([(1,4),(1,5),(9,13)])
        self.bplist22 = BasePairs([(1,4),(1,5),(1,13)])
        self.bplist222 = BasePairs([(1,13),(2,13),(3,13)])
        self.seq2 = 'GCACACACACACA'
        self.ref_bpseq2 = \
"""1 G 4
1 G 5
2 C 0
3 A 0
4 C 1
5 A 1
6 C 0
7 A 0
8 C 0
9 A 13
10 C 0
11 A 0
12 C 0
13 A 9
"""
        self.ref_bpseq22 = \
"""1 G 4
1 G 5
1 G 13
2 C 0
3 A 0
4 C 1
5 A 1
6 C 0
7 A 0
8 C 0
9 A 0
10 C 0
11 A 0
12 C 0
13 A 1
"""

        self.ref_bpseq222 = \
"""1 G 13
2 C 13
3 A 13
4 C 0
5 A 0
6 C 0
7 A 0
8 C 0
9 A 0
10 C 0
11 A 0
12 C 0
13 A 1
13 A 2
13 A 3
"""

        self.bplist3 = [(1,4),(1,5),(9,13)]
        self.seq3 = 'UCACACACACACA'

    def test_incorrect_init(self):
        """BpseqWriter.__init__ should not initalize nicely"""
        self.assertRaises(\
            BpseqWriterError, BpseqWriter, self.bplist3, self.seq3, self.path)

    def test_correct_init(self):
        """BpseqWriter.__init__ should initialize nicely"""
        bpseqwriter = BpseqWriter(self.bplist1, self.seq1, self.path)
        self.assertEqual(bpseqwriter.bplist, BasePairs([(1,4),(2,5),(9,13)]))
        self.assertEqual(bpseqwriter.seq, self.seq1)
        self.assertEqual(bpseqwriter.output_file.name, self.path)

    def test_calculateBpseq_bplist1(self):
        """BpseqWriter.calculateBpseq shoud calculate RNA ss of self.bpseq1"""
        bpseq_writer = BpseqWriter(self.bplist1, self.seq1, self.path)
        bpseq_writer.calculateBpseq()
        bpseq_writer.makeBpseqFile()
        self.assertEqual(bpseq_writer.bpseq, self.ref_bpseq1)
        self.assertEqual(parse_bp(open(self.path)).directed(), \
            list(self.bplist1.directed()))

    def test_calculateBpseq_bplist2(self):
        """BpseqWriter.calculateBpseq shoud calculate RNA ss of self.bpseq2"""
        bpseq_writer = BpseqWriter(self.bplist2, self.seq2, self.path)
        bpseq_writer.calculateBpseq()
        bpseq_writer.makeBpseqFile()
        self.assertEqual(bpseq_writer.bpseq, self.ref_bpseq2)
        #self.assertEqual(parse_bp(open(self.path)).directed(), \
            #list(self.bplist2.directed()))

    def test_calculateBpseq_bplist22(self):
        """BpseqWriter.calculateBpseq shoud calculate RNA ss of self.bpseq22"""
        bpseq_writer = BpseqWriter(self.bplist22, self.seq2, self.path)
        bpseq_writer.calculateBpseq()
        self.assertEqual(bpseq_writer.bpseq, self.ref_bpseq22)

    def test_calculateBpseq_bplist222(self):
        """BpseqWriter.calculateBpseq shoud calculate RNA ss of self.bpseq222"""
        bpseq_writer = BpseqWriter(self.bplist222, self.seq2, self.path)
        bpseq_writer.calculateBpseq()
        self.assertEqual(bpseq_writer.bpseq, self.ref_bpseq222)

    def test_makeBpseqFile(self):
        """BpseqWriter.makeBpseqFile should write self.bpseq to self.path"""
        bpseq_writer = BpseqWriter(self.bplist1, self.seq1, self.path)
        bpseq_writer.calculateBpseq()
        bpseq_writer.makeBpseqFile()
        self.assertEqual(open(self.path).read(), self.ref_bpseq1)

        bpseq_writer = BpseqWriter(self.bplist22, self.seq2, self.path)
        bpseq_writer.calculateBpseq()
        bpseq_writer.makeBpseqFile()
        self.assertEqual(open(self.path).read(), self.ref_bpseq22)

    def tearDown(self):
        """Cleans up leftovers"""
        os.close(self.fd)
        os.unlink(self.path)

if __name__ == '__main__':
    main()

