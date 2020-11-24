#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from unittest import TestCase, main
from wrappers import RASP


class RASPTests(TestCase):
    """ Tests for RASP"""
    def set_up(self):
        """Sets up test suite."""
        pass

    def test_potential_type_all_global_assessment(self):
        """RASP: Test all atom representation. Global assessment mode."""
        correct = -0.253625
        rasp = RASP('', '')
        try:
            path = os.path.dirname(__file__) \
                    and os.path.dirname(__file__) or "."
            result = rasp.run(path + os.sep + '1a9n.pdb', 'bbr', True)
        except Exception as e:
            print e
            result = ''
        finally:
            rasp.cleanup()
        self.assertEqual(correct, result)

    def test_potential_type_bbr_global_assessment(self):
        """RASP: Test bbr atom representation. Global assessment mode."""
        correct = -0.253625
        rasp = RASP('', '')
        try:
            path = os.path.dirname(__file__) \
                    and os.path.dirname(__file__) or "."
            result = rasp.run(path + os.sep + '1a9n.pdb', 'bbr', True)
        except Exception as e:
            print e
            result = ''
        finally:
            rasp.cleanup()
        self.assertEqual(correct, result)

    def test_potential_type_c3_global_assessment(self):
        """RASP: Test c3 atom representation. Global assessment mode."""
        correct = -0.0136507
        rasp = RASP('', '')
        try:
            path = os.path.dirname(__file__) \
                    and os.path.dirname(__file__) or "."
            result = rasp.run(path + os.sep + '1a9n.pdb', 'c3', True)
        except Exception as e:
            print e
            result = ''
        finally:
            rasp.cleanup()
        self.assertEqual(correct, result)

    def test_potential_type_all_profile(self):
        """RASP: Test all atom representation, profile mode."""
        correct = [(1, -775.03800000000001),
                   (2, -1164.22),
                   (3, -2054.1700000000001),
                   (4, -1601.1300000000001),
                   (5, -1619.45),
                   (6, -1177.2),
                   (7, -990.87400000000002),
                   (8, -781.63499999999999),
                   (9, -206.39699999999999),
                   (10, -147.00700000000001),
                   (11, -149.93700000000001),
                   (12, -142.125),
                   (13, -34.9236),
                   (14, -12.020300000000001),
                   (15, -376.87299999999999),
                   (16, -203.92699999999999),
                   (17, -1636.1500000000001),
                   (18, -930.91700000000003),
                   (19, -1194.3499999999999),
                   (20, -1363.79),
                   (21, -1542.6900000000001),
                   (22, -2767.9499999999998),
                   (23, -1288.26),
                   (24, -809.53800000000001)]
        rasp = RASP('', '')
        try:
            path = os.path.dirname(__file__) \
                    and os.path.dirname(__file__) or "."
            result = rasp.run( path + os.sep + '1a9n.pdb', 'all', False)
        except Exception as e:
            print e
            result = ''
        finally:
            rasp.cleanup()
        self.assertEqual(correct, result)

    def test_potential_type_bbr_profile(self):
        """RASP: Test bbr representation, profile mode."""
        correct = [(1, -295.09399999999999),
                   (2, -504.65899999999999),
                   (3, -819.01800000000003),
                   (4, -615.19100000000003),
                   (5, -654.56500000000005),
                   (6, -606.01199999999994),
                   (7, -464.71600000000001),
                   (8, -438.65899999999999),
                   (9, -284.89499999999998),
                   (10, -236.524),
                   (11, -256.18900000000002),
                   (12, -279.60599999999999),
                   (13, -237.35900000000001),
                   (14, -153.39599999999999),
                   (15, -296.69099999999997),
                   (16, -278.33600000000001),
                   (17, -774.43600000000004),
                   (18, -379.11700000000002),
                   (19, -624.29899999999998),
                   (20, -656.24900000000002),
                   (21, -585.47299999999996),
                   (22, -1231.79),
                   (23, -521.37300000000005),
                   (24, -341.39699999999999)]
        rasp = RASP('', '')
        try:
            path = os.path.dirname(__file__) \
                    and os.path.dirname(__file__) or "."
            result = rasp.run(path + os.sep + '1a9n.pdb', 'bbr', False)
        except:
            result = ''
        finally:
            rasp.cleanup()
        self.assertEqual(correct, result)

    def test_potential_type_c3_profile(self):
        """RASP: Test c3 representation, profile mode."""
        correct = [(1, -0.67069999999999996),
                   (2, -0.74909999999999999),
                   (3, -0.044400000000000002),
                   (4, -0.13270000000000001),
                   (5, -0.68279999999999996),
                   (6, 0.0247),
                   (7, 0.71060000000000001),
                   (8, -0.0088999999999999999),
                   (9, 3.4386100000000002),
                   (10, 0.1585),
                   (11, -0.30230000000000001),
                   (12, 0.31840000000000002),
                   (13, -0.075999999999999998),
                   (14, 0.15709999999999999),
                   (15, 0.071900000000000006),
                   (16, -0.049599999999999998),
                   (17, -0.1293),
                   (18, 0.0332),
                   (19, -0.75290000000000001),
                   (20, 217400000.0),
                   (21, 217400000.0),
                   (22, 217400000.0),
                   (23, -0.51900000000000002),
                   (24, -0.16750000000000001)]
        rasp = RASP('', '')
        try:
            path = os.path.dirname(__file__) \
                    and os.path.dirname(__file__) or "."
            result = rasp.run(path + os.sep + '1a9n.pdb', 'c3', False)
        except:
            result = ''
        finally:
            rasp.cleanup()
        self.assertEqual(correct, result)


if __name__ == '__main__':
    main()
