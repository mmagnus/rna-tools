#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test unit for SimRNA"""

from os import sep

from unittest import TestCase, main
from wrappers.SimRNA.SimRNA import SimRNA
from Config import WRAPPERS_PATH


class SimRNATests(TestCase):
    """ Tests for RASP"""
    def set_up(self):
        """SimRNA: Sets up test suite."""
        pass

    def test_0(self):
        """SimRNA: Only test"""
        correct = -162.732665
        simrna = SimRNA('', '')
        try:
            result = simrna.run(WRAPPERS_PATH + sep + 'SimRNA' + sep + 'test'
                                + sep + '1p5n_0_A-000001_corrected.pdb')
        except Exception as exception:
            print exception
            result = ''
        finally:
            simrna.cleanup()
        self.assertAlmostEqual(correct, result)

if __name__ == '__main__':
    main()
