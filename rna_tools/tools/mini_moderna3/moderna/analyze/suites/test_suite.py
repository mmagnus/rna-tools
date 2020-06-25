#!/usr/bin/env python
#
# test_suite.py
#
# Tests the suite module.
#

__author__ = "Kristian Rother, Raphael Bauer"
__credits__ = ["Marcin Domagalski","Magdalena Musielak", "Janusz Bujnicki", "Marie Curie"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

from unittest import main, TestCase
from suite import SuiteAngles, Suite
from suite_clusters import ClusterSet
from random import random
from numpy import array
from rna_tools.tools.mini_moderna3.moderna.PDB.PDBParser import PDBParser


class SuiteAngleTests(TestCase):
    def setUp(self):
        self.zero = SuiteAngles([0.0,0.0,0.0,0.0,0.0,0.0,0.0])
        self.a = SuiteAngles([1.0,1.0,1.0,1.0,1.0,1.0,1.0])
        self.b = SuiteAngles([1.0,2.0,3.0,4.0,5.0,6.0,7.0])
        self.c = SuiteAngles([5.0,5.0,3.0,4.0,5.0,6.0,5.0])
        self.d = SuiteAngles([1.5,1.0,1.0,1.5,1.0,1.0,1.0])
        self.e = SuiteAngles([1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0])
        
    def test_getitem(self):
        """indices 0 to 6 should return angles."""
        for i in range(7):
            self.assertEqual(self.zero[i],0.0)
            self.assertEqual(self.a[i],1.0)
            self.assertEqual(self.b[i],i+1.0)

    def test_set_angles(self):
        """changed angles should be visible for __getitem__"""
        self.zero.set_angles([3.2]*7)
        self.a.set_angles([0.1,2.1,4.1,6.1,8.1,10.1,12.1])
        for i in range(7):
            self.assertEqual(self.zero[i],3.2)
            self.assertEqual(self.a[i],i*2.0+0.1)
            
    def test_get_7d_diff(self):
        """should return indexable type with differences."""
        za = self.zero.get_7d_diff(self.a)
        #import pdb
        #pdb.set_trace()        
        az = self.a.get_7d_diff(self.zero)
        ZA_VALID = [1.0,1.0,1.0,1.0,1.0,1.0,1.0]
        bc = self.b.get_7d_diff(self.c)
        BC_VALID = [4.0,3.0,0.0,0.0,0.0,0.0,-2.0]        
        da = self.d.get_7d_diff(self.a)
        DA_VALID = [-0.5,0.0,0.0,-0.5,0.0,0.0,0.0]        
        for i in range(7):
            self.assertEqual(za[i],ZA_VALID[i])
            self.assertEqual(az[i],-ZA_VALID[i])
            self.assertEqual(bc[i],BC_VALID[i])
            self.assertEqual(da[i],DA_VALID[i])            

    def test_dotproduct(self):
        """Dotproducts for 4D and 7D should work"""
        za7 = self.zero.dotproduct(self.a, 7)
        self.assertEqual(za7, 0.0)
        aa7 = self.a.dotproduct(self.a, 7)
        self.assertEqual(aa7, 7.0)
        aa4 = self.a.dotproduct(self.a, 4)
        self.assertEqual(aa4, 4.0)        
        bc7 = self.b.dotproduct(self.c, 7)
        self.assertEqual(bc7, 136.0)        
        bc4 = self.b.dotproduct(self.c, 4)
        self.assertEqual(bc4, 60.0)        
        da7 = self.d.dotproduct(self.a, 7)
        self.assertEqual(da7, 8.0)        
        da4 = self.d.dotproduct(self.a, 4)
        self.assertEqual(da4, 4.5)
        ae7 = self.a.dotproduct(self.e, 7)
        self.assertEqual(ae7, 1.0)

    def test_get_hyperellipsoid_dist(self):
        weight_simple = [1.0]*7
        weight_complex = [0.1,0.1,0.1, 1.0,1.0,1.0, 10.0]
        zz7 = self.zero.get_hyperellipsoid_dist(self.zero,7,weight_simple)
        self.assertEqual(zz7, 0.0)
        zz4 = self.zero.get_hyperellipsoid_dist(self.zero,4,weight_simple)
        self.assertEqual(zz4, 0.0)
        za7 = self.zero.get_hyperellipsoid_dist(self.a,7,weight_simple)
        self.assertAlmostEqual(za7**3, 7.0)
        za4 = self.zero.get_hyperellipsoid_dist(self.a,4,weight_simple)
        self.assertAlmostEqual(za4**3, 4.0)
        aa7 = self.a.get_hyperellipsoid_dist(self.a,7,weight_complex)
        self.assertEqual(aa7, 0.0)
        za4 = self.zero.get_hyperellipsoid_dist(self.a,4,weight_complex)
        self.assertAlmostEqual(za4**3, 2002.0)
        za7 = self.zero.get_hyperellipsoid_dist(self.a,7,weight_complex)
        self.assertAlmostEqual(za7**3, 3003.001)
        bc4 = self.b.get_hyperellipsoid_dist(self.c,4,weight_simple)
        self.assertAlmostEqual(bc4**3, 27.0)
        bc7 = self.b.get_hyperellipsoid_dist(self.c,7,weight_simple)
        self.assertAlmostEqual(bc7**3, 99.0)
        bc4 = self.b.get_hyperellipsoid_dist(self.c,4,weight_complex)
        self.assertAlmostEqual(bc4**3, 27000.0)
        bc7 = self.b.get_hyperellipsoid_dist(self.c,7,weight_complex)
        self.assertAlmostEqual(bc7**3, 91000.008)
                    
                
class SuiteTests(TestCase):

    def test_assign_name_center(self):
        pass

if __name__== '__main__':
    main()

