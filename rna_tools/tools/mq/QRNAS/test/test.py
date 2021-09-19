#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
ABSPATH = os.path.abspath(__file__)

from unittest import TestCase, main
from wrappers.QRNA.QRNA import QRNA
from Config import WRAPPERS_PATH


class QRNATests(TestCase):
        """ Tests for RASP"""

        def setUp(self):
            """QRNA: Sets up test suite."""
            pass

        def test_electrostat(self):
            """QRNA: test with electrostatics
            """
            correct_energy = 29380.3125
            qrna = QRNA('GGGUGCUCAGUACGAGAGGAACCGCACCC', '1scl_0_A') ### TODO, you don't need a seq !!!
            try:
                energy = qrna.run(WRAPPERS_PATH + os.sep + 'QRNA' \
                        + os.sep + "test" + os.sep \
                        + "unmod_Val3_tRNA_model_si.pdb", electrostatics=True)
            except Exception as e:
                print e
                energy = ''
            finally:
                qrna.cleanup()
            self.assertAlmostEqual(correct_energy, energy)
            
        def test_noelectrostat(self):
            """QRNA: test without electrostatics
            """
            correct_energy = 29380.3125
            qrna = QRNA('GGGUGCUCAGUACGAGAGGAACCGCACCC', '1scl_0_A') ### TODO, you don't need a seq !!!
            try:
                energy = qrna.run(WRAPPERS_PATH + os.sep + 'QRNA' \
                        + os.sep + "test" + os.sep \
                        + "unmod_Val3_tRNA_model_si.pdb")
            except Exception as e:
                print e
                energy = ''
            finally:
                qrna.cleanup()
            self.assertAlmostEqual(correct_energy, energy)

if __name__ == '__main__':
        main()
