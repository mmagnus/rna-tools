#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test unit for lib.format"""

import os
import sys
import shutil
import uuid
from tempfile import tempdir
from unittest import TestCase, main

__file__ == 'test.py'

DIRNAME = os.path.dirname(__file__)
TEST_LINE = 'ATOM     14  C4\'   U A   5      25.964   2.534  -0.000  1.00  0.0               '

sys.path.append(os.path.abspath(os.path.join(DIRNAME, os.path.pardir)))
import SingleLineUtils
from PDBFile import PDBFile


class FormatlibTests(TestCase):
    """ Tests for formatlib"""
    def set_up(self):
        """Formatlib: Sets up test suite."""
        pass

    # SingleLineUtils tests
    def test_get_atom_code(self):
        result = SingleLineUtils.get_atom_code(TEST_LINE)
        self.assertEqual(result, 'C4\'')

    def test_get_atom_coords(self):
        result = SingleLineUtils.get_atom_coords(TEST_LINE)
        self.assertEqual(result, (5.964, 2.534, -0.0))

    def test_get_res_code(self):
        result = SingleLineUtils.get_res_code(TEST_LINE)
        self.assertEqual(result, '  U')

    def test_get_res_num(self):
        result = SingleLineUtils.get_res_num(TEST_LINE)
        self.assertEqual(result, 5)

    def test_set_line_bfactor(self):
        result = SingleLineUtils.set_line_bfactor(TEST_LINE, 10.5)
        correct = 'ATOM     14  C4\'   U A   5      25.964   2.534  -0.000  1.00 10.50              '
        self.assertEqual(result, correct)

    def test_set_atom_code(self):
        result = SingleLineUtils.set_atom_code(TEST_LINE, 'Xx')
        correct = 'ATOM     14  Xx    U A   5      25.964   2.534  -0.000  1.00  0.0               '
        self.assertEqual(result, correct)

    def test_set_res_code(self):
        result = SingleLineUtils.set_res_code(TEST_LINE, 'Xx')
        correct = 'ATOM     14  C4\'   XxA   5      25.964   2.534  -0.000  1.00  0.0               '
        self.assertEqual(result, correct)

    # tests for PDBFile
    def test_validate_pdb_true(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test0.in'))
        self.assertTrue(pdb_file.validate_pdb())

    def test_validate_pdb_false(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test.py'))
        self.assertFalse(pdb_file.validate_pdb())

    def test_detect_proteins(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'protein.pdb'))
        self.assertTrue(pdb_file.detect_proteins())

    def test_seq_from_pdb(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test0.in'))
        self.assertEqual(pdb_file.seq_from_pdb(), 'GAGCCGUAUGCGAUGAAAGUCGCACGUACGGUUC')

    def test_get_fasta(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test0.in'))
        self.assertEqual(pdb_file.get_fasta(name='test'), '>test\nGAGCCGUAUGCGAUGAAAGUCGCACGUACGGUUC\n')

    def test_remove_non_atoms(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test_nmr.in'))
        pdb_file.remove_non_atoms()
        with open(os.path.join(DIRNAME, 'test_nmr_atom.out')) as f:
            correct = f.read()
        self.assertEqual(pdb_file.pdb_string, correct)

    def test_resname_check_and_1to3(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test0.in'))
        pdb_file.resname_check_and_3to1()
        with open(os.path.join(DIRNAME, 'test0_3to1.out')) as f:
            correct = f.read()
        self.assertEqual(pdb_file.pdb_string, correct)

    def test_terminate_chains(self):
        # TODO: write it
        pass

    def test_remove_short_chains(self):
        # TODO: write it
        pass

    def test_check_and_add_P_at_start(self):
        # TODO: write it
        pass

    def test_set_residues_bfactor(self):
        # TODO: write it
        pass

    def test_count_models(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test_nmr.in'))
        self.assertEqual(pdb_file.count_models(), 10)

    def test_check_and_get_first_model(self):
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test_nmr.in'))
        pdb_file.check_and_get_first_model()
        with open(os.path.join(DIRNAME, 'test_nmr.out')) as f:
            correct = f.read()
        self.assertEqual(pdb_file.pdb_string, correct)

    # global tests
    def test_0(self):
        """Formatlib: first test"""
        tmpfile = os.path.join(tempdir or '/tmp/', str(uuid.uuid4()))
        shutil.copyfile(os.path.join(DIRNAME, 'test0.in'), tmpfile)
        pdb_file = PDBFile(pdb_path=tmpfile)
        pdb_file.pedantic_pdb()
        result = pdb_file.pdb_string
        os.remove(tmpfile)
        f = open(os.path.join(DIRNAME, 'test0.out'))
        correct = f.read()
        f.close()
        self.assertEqual(correct, result)

    def test_1(self):
        """Formatlib: second test"""
        pdb_file = PDBFile(pdb_path=os.path.join(DIRNAME, 'test1.in'))
        pdb_file.pedantic_pdb()
        result = pdb_file.pdb_string
        f = open(os.path.join(DIRNAME, 'test1.out'))
        correct = f.read()
        f.close()
        self.assertEqual(correct, result)


if __name__ == '__main__':
    main()
