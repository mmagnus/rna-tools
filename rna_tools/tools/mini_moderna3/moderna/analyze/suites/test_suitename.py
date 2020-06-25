#!/usr/bin/env python
#
# test_suitename.py
#
# Tests the suitename module.
#

__author__ = "Kristian Rother, Raphael Bauer"
__credits__ = ["Marcin Domagalski","Magdalena Musielak", "Janusz Bujnicki", "Marie Curie"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

from unittest import main, TestCase
from suitename import SuiteName, SuiteNameProcessor
from suite_clusters import ClusterSet
from random import random
from numpy import array
from Scientific.Geometry import Vector
from PDB.PDBParser import PDBParser
#from parse_dangle import *
import os,re

PATH = 'rnaDB2005/'

SUITE_EXAMPLES = (
('5j','ar0027H','B',17),
('7a','ar0041H','A',6),
('&a','pr0037Hb','B',163),
('1b','pr0113Hd','D',208),
('1[','pr0019Hb','B',658),
('7p','pr0033Hb','B',8),
('2o','pr0033Hb','B',5),
('6g','pr0122Hr','R',151),
('6j','pte003Hb','',975),
('1t','pte003Hb','B',907),
('5q','pte003Hb','B',973),
('8d','rr0009Hc','C',1062),
('6d','rr0082H09','',116),
('1m','rr0082H09','',1940),
('1L','rr0082H09','',1460),
('2h','rr0082H09','',2540),
('4n','rr0082H09','',767),
('0i','rr0082H09','',940),
('6n','rr0082H09','',1773),
('9a','rr0082H09','',2582),
('1g','rr0082H09','',1864),
('7d','rr0082H09','',636),
('3d','rr0082H09','',2118),
('2[','rr0082H09','',264),
('4b','rr0082H09','',247),
('0b','rr0082H09','',453),
('1o','rr0082H09','',1108),
('7r','rr0082H09','',262),
('2a','rr0082H09','',1711),
('4a','rr0082H09','',2485),
('0a','rr0082H09','',265),
('#a','rr0082H09','',1371),
('6p','rr0082H09','',1315),
('3b','rr0082H09','',904),
('1z','rr0082H09','',1771),
('4p','rr0096Hawx','',873),
('4d','tr0001H','',59),
('1f','tr0001H','',22),
('4g','ur0012Ha','A',226),
('4s','ur0026H','',2655),
('3a','urb016H','A',2),
('5d','ur0020H','A',9),
('1a','ur0020H','',11),
('1c','ur0020H','A',28),
('5z','ur0026H','',2654),
('1e','ur0035H','',2665),
)

class SuitenameTests(TestCase):

    def setUp(self):
        self.cs = ClusterSet()

    def get_original_seqs(self):
        """load suites calculated with dangle+suitename."""
        suiteseqs = {}
        for l in open('test_data/original_suites.txt'):
            name,seq = l.strip().split()
            suiteseqs[name] = seq
        return suiteseqs


    def test_assign_name_from_angles(self):
        suiteseqs = self.get_original_seqs()
        tested = 0
        errors = 0 
        for fn in os.listdir('test_data/dihedrals'):
            print fn
            seq = suiteseqs.get(fn[:-7],'no such structure')
            i = 0
            for angles in get_dihedrals('test_data/dihedrals/'+fn):
                s = SuiteName()
                s.set_angles(angles)
                try:
                    result = s.assign_name(self.cs,True,False,False,False)
                except:
                    result = '!!'
                shouldbe = seq[i:i+2]
                if shouldbe == '':
                    break # some sequences are too short, skip for this test.
                tested += 1
                if result!=shouldbe:
                    if not re.search('between-dom-sat',str(s)):
                        print seq[:i+4]
                        errors += 1
                        print result,'!=',shouldbe
                        print s
                i += 2
        print 'tested angles:',tested,'errors:',errors
        self.assertEqual(errors,0)

    def test_assign_name_center(self):
        """The centers of all clusters should match their names."""
        result = {}
        for bin in self.cs.get_bins():
            for cluster in self.cs.get_clusters_for_bin(bin):
                if cluster.wannabe: continue
                name = cluster.name
                s = SuiteName()
                angles = cluster.angles[:]
                angles[3] += 0.01
                s.set_angles(angles)
                result[name] = (s.assign_name(self.cs,True,False,False,False),s)
                if name != result[name][0]:
                        print "ASSIGNMENT FAILED: ",name, result[name]
                        print s
        for r in result:
            name, s = result[r]
            self.assertEqual(r, name)
            self.assertTrue(s.suiteness > 0.95)
            self.assertTrue(s.distance < 0.05)

    def get_examples(self):
        # pre-process the list to parse structures more efficiently
        s = {}
        for name,struc,chain,resi in SUITE_EXAMPLES:
            s.setdefault(struc,[])
            s[struc].append((name,chain,resi))
        # loop through structures
        for struc in s:
            parser = PDBParser()
            structure= parser.get_structure(struc,PATH+struc+'.pdb')
            # loop through each example
            for name,chain,resi in s[struc]:
                yield name, structure, chain, resi

            
    def itest_assign_name_examples(self):
        """The examples from the paper should all match."""
        result = {}
        for name, structure, chain_id, resi in self.get_examples():
            # pick the right chain
            if not chain_id:
                chain = structure[0].child_list[0]
            else:
                chain = None
                for c in structure[0].child_list:
                    if c.id == chain_id:
                        chain = c
            # get both residues
            resi1 = None
            resi2 = None
            for r in chain.child_list:
                if resi-1 == r.id[1]: resi1 = r
                elif resi == r.id[1]: resi2 = r
            self.assertTrue(resi1 and resi2)
            
            # now do the calculation
            s = SuiteName()
            s.set_residues(resi1,resi2)
            result[name] = s.assign_name(self.cs,True)
            if name != result[name]:
                print "\nASSIGNMENT FAILED: %s should be %s"%(result[name],name)
                print s
            else: print name,

        for r in result:
            self.assertEqual(result[r],r)

    def itest_get_string_for_pdbfile(self):
        """Suites should be same as those from original suitename program."""
        snp = SuiteNameProcessor()
        suiteseqs = self.get_original_seqs()
        
        # go through rna05 dataset
        pairs = []
        matches = 0
        for fn in os.listdir(PATH):
            seq = snp.get_string_for_pdbfile(PATH+fn)
            seq = re.sub('tt|oo','!!',seq)
            seq = re.sub('--','',seq)
            orig_seq = suiteseqs.get(fn[:-4],'no such structure')
            pairs.append((seq,orig_seq))
            if seq!=orig_seq:
                print '\nMISMATCH IN STRUCTURE',fn
                print seq
                print orig_seq
            else:
                matches += 1

        print '\nMATCHES %s/%s\n'%(matches,len(pairs))

        for seq,orig_seq in pairs:
            self.assertEqual(seq,orig_seq)

if __name__== '__main__':
    main()

