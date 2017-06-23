#!/usr/bin/env python
#-*-coding: utf-8-*-

from unittest import TestCase, main
from .secstruc import Secstruc, BasePairs, Partners, Vienna, solve_conflicts, \
ViennaStructure, ConflictInBasePairsError

class SecstrucTests(TestCase):

    def setUp(self):
        self.hairpin = Secstruc(HAIRPIN_SS)
        self.bulge = Secstruc(BULGE_SS)
        self.internal_loop = Secstruc(INTERNAL_LOOP_SS)
        self.junction = Secstruc(JUNCTION_SS)
        self.helix = Secstruc("((()))",[0,1,2,6,7,8])
        self.helix_shifted = Secstruc("((()))",[2,3,4,8,9,10])
        self.helix_a = Secstruc("(((",[0,1,2])
        self.helix_b = Secstruc(")))",[6,7,8])

    def test_eq(self):
        self.assertEqual(self.helix,Secstruc("((()))",[0,1,2,6,7,8]))
        self.assertNotEqual(self.helix, self.helix_a)
        self.assertNotEqual(self.helix_a, self.helix_b)
        self.assertNotEqual(self.helix, self.helix_shifted)
        self.assertEqual(self.hairpin, "((((...))))")

    def test_add(self):
        added = self.helix_a + self.helix_b
        self.assertEqual(added, self.helix)
        self.assertNotEqual(added, self.helix_shifted)

    def test_slice(self):
        a = self.helix[:3]
        self.assertEqual(a, self.helix_a)
        b = self.helix[3:]
        self.assertEqual(b, self.helix_b)
        self.assertEqual(self.helix[2:4], Secstruc("()",[2,6]))

    def test_getitem(self):
        self.assertEqual(self.helix[2], Secstruc("(",[2]))
        self.assertEqual(self.helix[3], Secstruc(")",[6]))
        
    def test_get_bp_array(self):
        s = Secstruc(".(())...((..)).((()))()..")
        self.assertEqual(s.get_bp_array(), [-1, 4, 3, 2, 1, -1, -1, -1, \
            13, 12, -1, -1, 9, 8, -1, 20, 19, 18, 17, 16, 15, 22, 21, -1, -1])
        
    def test_get_type(self):
        """Should recognize substructures."""
        self.assertEqual(Secstruc('()').get_type(),None)
        self.assertEqual(Secstruc('((()))').get_type(),'helix')
        self.assertEqual(Secstruc('((())))').get_type(),None)
        self.assertEqual(Secstruc('(())').get_type(),'helix')        
        self.assertEqual(Secstruc('...(').get_type(),"5'-overhang")
        self.assertEqual(Secstruc(').').get_type(),"3'-overhang")
        self.assertEqual(Secstruc('(....)').get_type(),'loop')
        self.assertEqual(Secstruc('(().)').get_type(),'bulge')
        self.assertEqual(Secstruc('(.())').get_type(),'bulge')
        self.assertEqual(Secstruc('(..())').get_type(),'bulge')
        self.assertEqual(Secstruc('(..().())').get_type(),'junction')
        self.assertEqual(Secstruc('(()()())').get_type(),'junction')

    def test_find_overhang5(self):
        """Should part the object into an overhang secstruc."""
        # 5'-overhangs
        overhang = Secstruc("..((()))").find_overhang5()
        self.assertEqual(overhang, Secstruc('..(',[0,1,2]))

    def test_find_overhang3(self):
        """Should part the object into an overhang secstruc."""
        # 3'-overhangs
        overhang = Secstruc("(((..)))....").find_overhang3()
        self.assertEqual(overhang, Secstruc(')....',[7,8,9,10,11]))
        
    def test_find_helices(self):
        """Should return all helices from a secstruc."""
        # simple helix
        helices = Secstruc("...(())..").find_helices()
        self.assertEqual(len(helices),1)
        self.assertEqual(helices[0], Secstruc('(())', [3,4,5,6]))
        
        # multiple helices
        helices = Secstruc(".(())...((..)).((()))()..").find_helices()
        self.assertEqual(len(helices),4)
        
        # full helix in one.
        helices = Secstruc("((()))").find_helices()
        self.assertEqual(helices[0], Secstruc('((()))',[0,1,2,3,4,5]))
        
        # two helices with bulge in between
        helices = Secstruc("((((((..))).)))").find_helices()
        self.assertEqual(len(helices),2)
        self.assertEqual(helices[0], Secstruc('((()))',[0,1,2,12,13,14]))
        self.assertEqual(helices[1], Secstruc('((()))',[3,4,5,8,9,10]))

    def test_find_loops(self):
        """Should return all loops from a secstruc."""
        loops = Secstruc("(((..))).").find_loops()
        self.assertEqual(loops[0], Secstruc('(..)',[2,3,4,5]))

    def test_find_junctions(self):
        """Should return all junctions from a secstruc."""
        junctions = Secstruc('(()()())').find_junctions()
        self.assertEqual(junctions[0], Secstruc('(()()())'))
        # helices should be no junctions
        junctions = Secstruc('((()))').find_junctions()
        self.assertEqual(len(junctions), 0)

    def test_find_bulges(self):
        """Should return all bulges from a secstruc."""
        bulges = Secstruc(BULGE_SS).find_bulges()
        self.assertEqual(len(bulges),1)
        self.assertEqual(bulges[0], Secstruc('(().)',[3,4,11,12,13]))
        # also check bulge in other direction
        bulges = Secstruc("((.((..))))").find_bulges()
        self.assertEqual(bulges[0], Secstruc('(.())',[1,2,3,8,9]))
        
    def test_find_internal_loop(self):
        """Internal loops should be recognized as bulges."""
        loops = Secstruc(INTERNAL_LOOP_SS).find_bulges()
        self.assertEqual(len(loops),1)
        self.assertEqual(loops[0], Secstruc('(....()..)',\
        [2,3,4, 5, 6, 7, 15, 16, 17, 18]))
        
    def test_find_secstruc_elements(self):
        """Should divide a sec struct into substructures and their indices."""
        result = Secstruc(HAIRPIN_SS).find_secstruc_elements()
        self.assertTrue(Secstruc('(...)',[3,4,5,6,7]) in result)
        self.assertTrue(Secstruc('(((())))',[0,1,2,3,7,8,9,10]) in result)
        self.assertEqual(len(result), 2)

    def test_find_secstruc_elements_bulge(self):
        result = Secstruc(BULGE_SS).find_secstruc_elements()
        self.assertTrue(Secstruc('(....)',[5,6,7,8,9,10]) in result)
        self.assertTrue(Secstruc('(().)',[3,4,11,12,13]) in result)
        self.assertTrue(Secstruc('(())',[4,5,10,11]) in result)
        self.assertTrue(Secstruc('(((())))',[0,1,2,3,13,14,15,16]) in result)
        self.assertEqual(len(result), 4)

    def test_find_secstruc_elements_junction(self):
        result = Secstruc(JUNCTION_SS).find_secstruc_elements()
        self.assertTrue(Secstruc('..(',[0,1,2]) in result)
        self.assertTrue(Secstruc('((()))',[2,3,4,27,28,29]) in result)
        self.assertTrue\
            (Secstruc('(..().())',[4,5,6,7,17,18,19,26,27]) in result)
        self.assertTrue(Secstruc('(((())))',[7,8,9,10,14,15,16,17]) in result)
        self.assertTrue(Secstruc('(...)',[10,11,12,13,14]) in result)
        self.assertTrue(Secstruc('(())',[19,20,25,26]) in result)
        self.assertTrue(Secstruc('(....)',[20,21,22,23,24,25]) in result)
        self.assertEqual(len(result), 7)
        
HAIRPIN_SEQ = "AAAACCCUUUU"
HAIRPIN_SS  = "((((...))))"

BULGE_SEQ = "GGCCGGAAAACCUGGCC"
BULGE_SS  = "((((((....)).))))"

INTERNAL_LOOP_SEQ = "GGGAAAAUUUCCCAAAGGCCC"
INTERNAL_LOOP_SS = "(((....(((...)))..)))"

JUNCTION_SEQ = "AAUCGAAUAUAGCCUAUAGCCAUGCGGCGA"
JUNCTION_SS  = "..(((..((((...)))).((....)))))"
###############################################################################
##                                                                           ##
## CLASSES ViennaStructureTests, RnaAlphabet, Rna WERE COPIED                ##
## FROM rna2d.py MODULE FROM PyCogent WITHOUT ANY MAJOR CHANGES.             ##
## THE ONLY DIFFERENCE FROM PYCOGENT IS ViennaStructureTests WHICH DOES NOT  ##
## HAVE self.test_toTree() TEST.                                             ##
## (PyCogent ver. 1.3.1)                                                     ##
##                                                                           ##
###############################################################################

class PartnersTests(TestCase):
    """Tests for Partners object"""

    def setUp(self):
        self.empty = Partners([None] * 6)

    def test_init(self):
        """Partners should init with empty list and stay free of conflicts
        """
        self.assertEqual(Partners([]),[])
        self.assertEqual(self.empty, [None, None, None, None, None, None])
        self.assertRaises(ValueError, self.empty.__setitem__, 2, (2,))
        
    def test_set_1(self):
        """Setting partners should work, part 1
        """
        self.empty[2] = (3,)
        self.assertEqual(self.empty, [None, None, (3,), (2,), None, None])
        
    def test_set_2(self):
        """Setting partners should work, part 2
        """        
        self.empty[3] = (4, 5)
        self.assertEqual(self.empty,[None, None, None, (4, 5), (3,), (3,)])
        
    def test_set_3(self):
        """Setting partners should work, part 3
        """
        self.empty[3] = (5,)
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])
    
        self.empty[1] = None
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])
        
        self.empty[1] = (2,)
        self.assertEqual(self.empty, [None, (2,), (1,), (5,), None, (3,)])
        
    def test_set_4(self):
        """Setting partners should work, part 4
        """
        self.empty[3] = (5,)
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])
    
        self.empty[1] = None
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])
        
        self.empty[1] = (2,)
        self.assertEqual(self.empty, [None, (2,), (1,), (5,), None, (3,)])
        
        self.empty[1] = None
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])
        
    def test_set_5(self):
        """Setting partners should work, part 5
        """
        self.empty[3] = (5,)
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])
        
        self.empty[3] = (4,) 
        self.assertEqual(self.empty, [None, None, None, (4,), (3,), None])

        self.empty[3] = (5,)
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])

        self.empty[3] = (4, 5)
        self.assertEqual(self.empty, [None, None, None, (4, 5), (3,), (3,)])
        
    def test_set_6(self):
        """Setting partners should work, part 6
        """
        self.empty[3] = (5,)
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])

        self.empty[2] = (4,)
        self.assertEqual(self.empty, [None, None, (4,), (5,), (2,), (3,)])
        
    def test_set_7(self):
        """Setting partners should work, part 7
        """
        self.empty[3] = (5,)
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])

        self.empty[2] = (4,5)
        self.assertEqual(self.empty, [None, None, (4,5), (5,), (2,), (2,3)])
        
        self.empty[2] = (0,4,5)
        self.assertEqual(self.empty, [(2,), None, (0,4,5), (5,), (2,), (2,3)])

        self.empty[2] = None
        self.assertEqual(self.empty, [None, None, None, (5,), None, (3,)])
        
    def test_toPairs(self):
        """Partners toPairs() should return a BasePairs object
        """
        partners = Partners([None, None, None, (5,), None, (3,)])
        bps = partners.toPairs()
        ref_bps = BasePairs([(3,5)])
        self.assertEqual(bps, ref_bps)
        
        partners = Partners([None, None, (4,5), (5,), (2,), (2,3)])
        bps = partners.toPairs()
        ref_bps = BasePairs([(2,4), (2,5), (3,5)])
        self.assertEqual(bps, ref_bps)

        partners = Partners([(2,), None, (0,4,5), (5,), (2,), (2,3)])
        bps = partners.toPairs()
        ref_bps = BasePairs([(0,2), (2,4), (2,5), (3,5)])
        self.assertEqual(bps, ref_bps)
        
        partners = Partners([None, None, (4,), (5,), (2,), (3,)])
        bps = partners.toPairs()
        ref_bps = BasePairs([(2,4), (3,5)])
        self.assertEqual(bps, ref_bps)
        
        partners = Partners([None, (2,), (1,), (5,), None, (3,)])
        bps = partners.toPairs()
        ref_bps = BasePairs([(1,2), (3,5)])
        self.assertEqual(bps, ref_bps)
        
        partners = Partners([None, None, None, (4, 5), (3,), (3,)])
        bps = partners.toPairs()
        ref_bps = BasePairs([(3,4), (3,5)])
        self.assertEqual(bps, ref_bps)

        partners = Partners([None, None, None])
        bps = partners.toPairs()
        ref_bps = BasePairs([])
        self.assertEqual(bps, ref_bps)
        
class ViennaStructureTests(TestCase):
    """Test that ViennaStructure methods and properties work.
    
    These tests should work for every StructureString.
    """
    def setUp(self):
        """SetUp for all ViennaStructureTests"""
        self.Energy1 = 0.0
        self.Energy2 = -1e-02
        self.NoPairsStr = '.....'
        self.NoPairs = ViennaStructure('.....',self.Energy1)
        self.OneHelixStr = '((((()))))'
        self.OneHelix = ViennaStructure('((((()))))',self.Energy2)
        self.TwoHelixStr = '((.))(()).'
        self.TwoHelix = ViennaStructure('((.))(()).',self.Energy1)
        self.ThreeHelixStr = '(((((..))..(((..)))..)))'
        self.ThreeHelix = ViennaStructure('(((((..))..(((..)))..)))')
        self.EmptyStr = ''
        self.ManyHelicesStr = '(..(((...)).((.(((((..))).)))..((((..))))))...)'
        self.EndsStr = '..(.)..'   #has trailing bases at ends
        self.InternalStr = '(((...)))..(()).' #has internal non-nested region
        self.Pseudoknots = ViennaStructure('(.((..(.(..ab.cde..(.EDC.BA(.'+\
                                           '(((()))).).(...(....)...).fghi'+\
                                           '..)....IHGF).)..)).)...')
        self.pseudoknots_bps = BasePairs([(1, 79), (3, 77), (4, 76), (7, 73),
            (9, 71), (12, 27), (13, 26), (15, 24), (16, 23), (17, 22), (20, 62),
            (28, 39), (30, 37), (31, 36), (32, 35), (33, 34), (41, 54),
            (45, 50), (56, 70), (57, 69), (58, 68), (59, 67)])

    def test_len(self):
        """ViennaStructure len() should match structure length"""
        self.assertEqual(len(self.TwoHelix), 10)
        #Do you want the init possible with Vienna() for empty str???
        self.assertEqual(len(Vienna('')), 0)

    def test_getitem(self):
        """ViennaStructure struct[index] should return char at index in struct"""
        self.assertEqual(self.NoPairs[0], self.NoPairsStr[0])  #dot
        self.assertEqual(self.OneHelix[0], self.OneHelixStr[0]) #close pair
        #negative index, end
        self.assertEqual(self.OneHelix[-1], self.OneHelixStr[-1])
        #middle of sequence
        self.assertEqual(self.TwoHelix[4], self.TwoHelixStr[4]) 
                
    def test_getslice(self):
        """ViennaStructure struct[a:b] should return slice from a to b"""
        self.assertEqual(self.OneHelix[0:3], self.OneHelixStr[0:3])
        self.assertEqual(self.OneHelix[0:0], self.OneHelixStr[0:0])
        self.assertEqual(self.OneHelix[:], self.OneHelixStr[:])
        self.assertEqual(self.OneHelix[4:7], self.OneHelixStr[4:7])
        self.assertEqual(self.OneHelix[5:], self.OneHelixStr[5:])
        
    def test_str(self):
        """ViennaStructure str() should print structure and energy, if known"""
        self.assertEqual(str(self.NoPairs),  '..... (0.0)')
        self.assertEqual(str(self.OneHelix), '((((())))) (-0.01)')
        self.assertEqual(str(self.TwoHelix), '((.))(()). (0.0)')
        self.assertEqual(str(self.ThreeHelix), '(((((..))..(((..)))..)))')
        
    def test_toPartners(self):
        """ViennaStructure toPartners() should return Partners object"""
        self.assertEqual(self.NoPairs.toPartners(), [None]*5)
        self.assertEqual(self.OneHelix.toPartners(),[(9,), (8,), (7,), (6,),
            (5,), (4,), (3,), (2,), (1,), (0,)])
        self.assertEqual(self.TwoHelix.toPartners(),
                         [(4,), (3,), None, (1,), (0,), (8,), (7,),
                            (6,), (5,), None])
        self.assertEqual(self.ThreeHelix.toPartners(),
                         [(23,), (22,), (21,), (8,), (7,), None, None, (4,),
                            (3,), None, None, (18,), (17,), (16,), None, None,
                            (13,), (12,), (11,), None, None, (2,), (1,), (0,)])

    def test_toPairs(self):
        """ViennaStructure toPairs() should return Pairs object"""
        self.assertEqual(self.NoPairs.toPairs(), [])
        for item in self.OneHelix.toPairs():
            self.assertTrue(item in  [(1,10),(2,9),(3,8),(4,7),(5,6)])
        for item in self.TwoHelix.toPairs():
            self.assertTrue(item in [(1,5),(2,4),(6,9),(7,8)])
        for item in self.ThreeHelix.toPairs():
            self.assertTrue(item in \
            [(1,24),(2,23),(3,22),(4,9),(5,8),(12,19),(13,18),(14,17)])
        self.assertEqual(self.Pseudoknots.toPairs(), [(1, 79), (3, 77), \
            (4, 76), (7, 73), (9, 71), (12, 27), (13, 26), (15, 24), (16,23),\
            (17, 22), (20, 62), (28, 39), (30, 37), (31, 36), (32, 35), \
            (33, 34), (41, 54), (45, 50), (56, 70), (57, 69), (58, 68), \
            (59, 67)])
        
    def test_toPairs_vienna_conversion(self):
        """toPairs.toVienna() should generate the same vienna struct
        """
        self.assertEqual(self.OneHelix.toPairs().toVienna(len(self.OneHelix)),
                         self.OneHelix)
        self.assertEqual(self.TwoHelix.toPairs().toVienna(len(self.TwoHelix)),
                         self.TwoHelix)
        self.assertEqual(\
            self.ThreeHelix.toPairs().toVienna(len(self.ThreeHelix)),
            self.ThreeHelix)
        self.assertEqual(ViennaStructure(self.ManyHelicesStr).toPairs().\
                         toVienna(len(self.ManyHelicesStr)),
                         self.ManyHelicesStr)

    def test_toPairs_pseudoknotted_vienna_conversion(self):
        """toPairs.toVienna() should generate the same pseudoknotted vienna
        """
        self.assertEqual(self.Pseudoknots.toPairs(), self.pseudoknots_bps)

        self.assertEqual(\
            self.pseudoknots_bps.toVienna(len(self.Pseudoknots)),
            '(.((..(.(..((.(((..[.))).))(.(((()))).).(...(....).'+\
            '..).{{{{..]....}}}}).)..)).)...') # the same secondary structure
                                               # as self.Pseudoknots, the only
                                               # difference is notation, but
                                               # logic is the same !!!

        self.assertEqual(self.Pseudoknots.toPairs().\
                         toVienna(len(self.Pseudoknots)),
                         '(.((..(.(..((.(((..[.))).))(.(((()))).).(...(....).'+\
                         '..).{{{{..]....}}}}).)..)).)...') # the same secondary
                                                            # structure as
                                                            # self.Pseudoknots,
                                                            # the only
                                                            # difference is
                                                            # notation, but
                                                            # logic is the
                                                            # same !!!

class RnaAlphabet(object):
    Pairs = {
    ('A','U'): True,    #True vs False for 'always' vs 'sometimes' pairing
    ('C','G'): True,
    ('G','C'): True,
    ('U','A'): True,
    ('G','U'): False,
    ('U','G'): False,
}

class Rna(str):
    Alphabet = RnaAlphabet
###############################################################################
##                                                                           ##
###############################################################################

class BasePairsTests(TestCase):
    """Tests for BasePairs class from Secstruc module"""
    
    def test_identifyNestedBasepairs(self):
        """BasePairs.identifyNestedBasepairs identify nested base pairs"""
        self.assertEqual(self.no_pseudoknots.identifyNestedBasepairs(),\
        {BasePairs([(1,12)]) : [(2,9)], BasePairs([(13,15)]): []})
        
        self.assertEqual(self.two_pseudoknots_1.identifyNestedBasepairs(),\
        {BasePairs([(1,12)]) : [], BasePairs([(11,13)]) : [], \
         BasePairs([(14,16)]) : [], BasePairs([(15,17)]) : []})
        
        self.assertEqual(self.two_pseudoknots_2.identifyNestedBasepairs(),\
        {BasePairs([(1,12)]) : [], BasePairs([(9,13)]) : [], \
         BasePairs([(10,14)]) : [], BasePairs([(15,17)]) : []})
        
        self.assertEqual(self.two_pseudoknots_3.identifyNestedBasepairs(),\
        {BasePairs([(1,12)]) : [], BasePairs([(9,13)]) : [], \
         BasePairs([(10,19)]) : [], BasePairs([(20,117)]) : []})

        self.assertEqual(self.three_pseudoknots.identifyNestedBasepairs(),\
        {BasePairs([(1,12)]) : [(6,10)], \
         BasePairs([(9,20)]) : [(11,13)], \
         BasePairs([(14,16)]) : [], \
         BasePairs([(15,17)]): []})
        
        self.assertEqual(\
            self.bplist_with_one_long_pseudoknot.identifyNestedBasepairs(),\
            {BasePairs([(2,24)]) : [(3,23), (4,22), (5,21), (6,20), (7,19), \
             (8,18)], BasePairs([(11,45)]) : [(12,44), (13,43), (14,42), \
             (15,41), (16,40), (17,39)]})
            
        self.assertEqual(self.bplist_with_one_gapped_pseudoknot.\
            identifyNestedBasepairs(),\
            {BasePairs([(2,24)]) : [(3,23), (4,22), (5,21), (6,20), (7,19), \
             (8,18)], BasePairs([(11,45)]) : [(12,44), (13,43), \
             (15,41), (16,40), (17,39)]})
            
        self.assertEqual(self.bplist_with_two_long_pseudoknots_1.\
            identifyNestedBasepairs(),\
            {BasePairs([(2,22)]) : [(3,21), (4,20), (5,19), (6,18), (7,17), \
             (8,16)], BasePairs([(11,26)]) : [(12,25), (13,24)], \
             BasePairs([(29,41)]): [(30,40), (31,39)], \
             BasePairs([(34,46)]) : [(35,45), (36,44)]})
            
        self.assertEqual(self.bplist_with_two_long_pseudoknots_2.\
            identifyNestedBasepairs(),\
            {BasePairs([(2,22)]) : [(3,21), (4,20), (5,19), (6,18), (7,17), \
             (8,16)], BasePairs([(11,46)]) : [(12,45), (13,44), (24,36),\
             (25,35), (26,34)], BasePairs([(29,41)]): [(30,40), (31,39)]})
            
        self.assertEqual(self.bplist_with_two_long_pseudoknots_3.\
            identifyNestedBasepairs(),\
            {BasePairs([(2,26)]) : [(3,25), (4,24), (5,23), (6,22), (7,21), \
             (8,20)], BasePairs([(11,37)]) : [(12,36), (13,35)], \
             BasePairs([(30,41)]): [(31,40), (32,39)]})
            
        self.assertEqual(self.bplist_with_two_gapped_pseudoknots_1.\
            identifyNestedBasepairs(),\
            {BasePairs([(2,22)]) : [(3,21), (4,20), (5,19), (6,18), (7,17), \
             (8,16)], BasePairs([(11,46)]) : [(13,44), (24,36), (25, 35), \
             (26,34)], BasePairs([(29,41)]): [(30,40), (31,39)]})

        self.assertEqual(self.bplist_with_two_gapped_pseudoknots_2.\
            identifyNestedBasepairs(),\
            {BasePairs([(2,22)]) : [(3,21), (4,20), (6,18), (7,17), \
             (8,16)], BasePairs([(11,46)]) : [(12,45), (13,44), (24,36), \
             (26,34)], BasePairs([(29,41)]): [(30,40), (31,39)]})
            
        self.assertEqual(self.bplist_with_conflicts.identifyNestedBasepairs(),\
            {BasePairs([(1, 77)]): [(2, 72), (3, 71), (4, 70),(5, 69),(6, 68),\
             (10, 48),(11, 47),(12, 46),(17, 42),(18, 41),(19, 40),(20, 39),\
             (21, 38),(22, 37),(23, 35),(24, 34),(25, 33),(26, 32)],\
             BasePairs([(50, 66)]): [(51, 65), (52, 64), (53, 63), (54, 62)]})
            
    def test_hasPseudoknots(self):
        """BasePairs.hasPseudoknots should return True if pseudoknot exists"""
        self.assertFalse(self.no_pseudoknots.hasPseudoknots())
        self.assertFalse(self.bplist_with_conflicts.hasPseudoknots())
        self.assertTrue(self.one_pseudoknot.hasPseudoknots())
        self.assertTrue(self.two_pseudoknots_1.hasPseudoknots())
        self.assertTrue(self.two_pseudoknots_2.hasPseudoknots())
        self.assertTrue(self.two_pseudoknots_3.hasPseudoknots())
        self.assertTrue(self.three_pseudoknots.hasPseudoknots())
        self.assertTrue(self.bplist_with_one_long_pseudoknot.hasPseudoknots())
        self.assertTrue(\
            self.bplist_with_one_gapped_pseudoknot.hasPseudoknots())
        self.assertTrue(\
            self.bplist_with_two_long_pseudoknots_1.hasPseudoknots())
        self.assertTrue(\
            self.bplist_with_two_long_pseudoknots_2.hasPseudoknots())
        self.assertTrue(\
            self.bplist_with_two_long_pseudoknots_3.hasPseudoknots())
        self.assertTrue(\
            self.bplist_with_two_gapped_pseudoknots_1.hasPseudoknots())
        self.assertTrue(\
            self.bplist_with_two_gapped_pseudoknots_2.hasPseudoknots())
            
    def test_getPseudoknots(self):
        """BasePairs.getPseudoknots should correctly get pseudoknots"""
        self.assertEqual(self.no_pseudoknots.getPseudoknots(), [])
        self.assertEqual(self.bplist_with_conflicts.getPseudoknots(), [])
        self.assertEqual(self.one_pseudoknot.getPseudoknots(), [(11, 13)])
        self.assertEqual(self.two_pseudoknots_1.getPseudoknots(),
                         [(11, 13), (15, 17)])
        self.assertEqual(self.two_pseudoknots_2.getPseudoknots(),
                         [(9, 13), (10, 14)])
        self.assertEqual(self.two_pseudoknots_3.getPseudoknots(),
                         [(9, 13), (10, 19)])
        self.assertEqual(self.three_pseudoknots.getPseudoknots(),
                         [(9, 20), (11, 13), (15, 17)])
        self.assertEqual(\
            self.bplist_with_one_long_pseudoknot.getPseudoknots(),
            [(11, 45), (12, 44), (13, 43), (14, 42), (15, 41), (16, 40),
             (17, 39)])
        self.assertEqual(\
            self.bplist_with_one_gapped_pseudoknot.getPseudoknots(),
            [(11, 45), (12, 44), (13, 43), (15, 41), (16, 40), (17, 39)])
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_1.getPseudoknots(),
            [(11, 26), (12, 25), (13, 24), (34, 46), (35, 45), (36, 44)])
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_2.getPseudoknots(),
            [(11, 46), (12, 45), (13, 44), (29, 41), (30, 40), (31, 39)])
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_3.getPseudoknots(),
            [(11, 37), (12, 36), (13, 35), (30, 41), (31, 40), (32, 39)])
            #{(8, 20): [(11, 37), (12, 36), (13, 35)]}) # CHECK!
        self.assertEqual(\
            self.bplist_with_two_gapped_pseudoknots_1.getPseudoknots(),
            [(11, 46), (13, 44), (29, 41), (30, 40), (31, 39)])
        self.assertEqual(\
            self.bplist_with_two_gapped_pseudoknots_2.getPseudoknots(),
            [(11, 46), (12, 45), (13, 44), (29, 41), (30, 40), (31, 39)])
        
        self.assertEqual(BasePairs([(1, 79), (3, 77), (4, 76), (7, 73),
                                    (9, 71), (12, 27), (13, 26), (15, 24),
                                    (16, 23), (17, 22), (20, 62), (28, 39),
                                    (30, 37), (31, 36), (32, 35), (33, 34),
                                    (41, 54), (45, 50), (56, 70), (57, 69),
                                    (58, 68), (59, 67)]).getPseudoknots(),
                         [(20, 62), (56, 70), (57, 69), (58, 68), (59, 67)])
        
    def test_pseudoknots(self):
        """BasePairs.pseudoknots should correctly get pseudoknots"""
        self.assertEqual(self.no_pseudoknots.pseudoknots, [])
        self.assertEqual(self.bplist_with_conflicts.pseudoknots, [])
        self.assertEqual(self.one_pseudoknot.pseudoknots, [(11, 13)])
        self.assertEqual(self.two_pseudoknots_1.pseudoknots,
                         [(11, 13), (15, 17)])
        self.assertEqual(self.two_pseudoknots_2.pseudoknots,
                         [(9, 13), (10, 14)])
        self.assertEqual(self.two_pseudoknots_3.pseudoknots,
                         [(9, 13), (10, 19)])
        self.assertEqual(self.three_pseudoknots.pseudoknots,
                         [(9, 20), (11, 13), (15, 17)])
        self.assertEqual(\
            self.bplist_with_one_long_pseudoknot.pseudoknots,
            [(11, 45), (12, 44), (13, 43), (14, 42), (15, 41), (16, 40),
             (17, 39)])
        self.assertEqual(\
            self.bplist_with_one_gapped_pseudoknot.pseudoknots,
            [(11, 45), (12, 44), (13, 43), (15, 41), (16, 40), (17, 39)])
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_1.pseudoknots,
            [(11, 26), (12, 25), (13, 24), (34, 46), (35, 45), (36, 44)])
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_2.pseudoknots,
            [(11, 46), (12, 45), (13, 44), (29, 41), (30, 40), (31, 39)])
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_3.pseudoknots,
            [(11, 37), (12, 36), (13, 35), (30, 41), (31, 40), (32, 39)])
            #{(8, 20): [(11, 37), (12, 36), (13, 35)]}) # CHECK!
        self.assertEqual(\
            self.bplist_with_two_gapped_pseudoknots_1.pseudoknots,
            [(11, 46), (13, 44), (29, 41), (30, 40), (31, 39)])
        self.assertEqual(\
            self.bplist_with_two_gapped_pseudoknots_2.pseudoknots,
            [(11, 46), (12, 45), (13, 44), (29, 41), (30, 40), (31, 39)])
        
        self.assertEqual(BasePairs([(1, 79), (3, 77), (4, 76), (7, 73),
                                    (9, 71), (12, 27), (13, 26), (15, 24),
                                    (16, 23), (17, 22), (20, 62), (28, 39),
                                    (30, 37), (31, 36), (32, 35), (33, 34),
                                    (41, 54), (45, 50), (56, 70), (57, 69),
                                    (58, 68), (59, 67)]).pseudoknots,
                         [(20, 62), (56, 70), (57, 69), (58, 68), (59, 67)])
            
    def test_toVienna_wrong_length(self):
        """BasePairs.toVienna should raise IndexError for wrong length"""
        self.assertRaises(IndexError, self.no_pseudoknots.toVienna, 10, -1)
        self.assertRaises(IndexError, self.no_pseudoknots.toVienna, 0, 0)
    
    def test_toVienna_wrong_offset(self):
        """BasePairs.toVienna should raise IndexError for wrong offset"""
        self.assertRaises(IndexError, self.one_pseudoknot.toVienna, 15, -8)
        self.assertRaises(IndexError, self.one_pseudoknot.toVienna, 20, 20)

    def test_toPartners(self):
        """BasePairs.toPartners should return a Partners object"""
        a = BasePairs([(1,5),(3,4),(6,9),(7,8)]) #normal
        b = BasePairs([(0,4),(2,6)]) #pseudoknot
        c = BasePairs([(1,6),(3,6),(4,5)]) #conflict

        self.assertEqual(a.toPartners(10),
                         [None, (5,), None, (4,), (3,), (1,), (9,), (8,), (7,),
                          (6,)])
        self.assertEqual(a.toPartners(13,3),\
        [None, None, None, None, (8,), None, (7,), (6,), (4,), (12,), (11,),
         (10,), (9,)])
        assert isinstance(a.toPartners(10),Partners)
        self.assertEqual(b.toPartners(7),
                         [(4,), None, (6,), None, (0,), None, (2,)])
        #self.assertRaises(ValueError,c.toPartners,7, strict=True)
        self.assertEqual(c.toPartners(7),
                         [None, (6,), None, (6,), (5,), (4,), (1,3)])

        #raises an error when try to insert something at non-existing indices
        self.assertRaises(IndexError, c.toPartners, 0)

    def test_toVienna_OK(self):
        """BasePairs.toVienna should produce Vienna seconadry structure"""
        self.assertEqual(self.no_pseudoknots.toVienna(16, -1),
                         '((......)..)(.).')
        self.assertEqual(self.no_pseudoknots.toVienna(16, -1).toPairs(),
                         self.no_pseudoknots)
        
        self.assertEqual(self.one_pseudoknot.toVienna(16, -1),
                         '(.........[)](.)')
        self.assertEqual(self.one_pseudoknot.toVienna(16, -1).toPairs(),
                         self.one_pseudoknot)
        
        self.assertEqual(self.two_pseudoknots_1.toVienna(18, -1),
                         '(.........[)]([)].')
                         
        self.assertEqual(self.two_pseudoknots_1.toVienna(18, -1).toPairs(),
                         self.two_pseudoknots_1)

        self.assertEqual(self.two_pseudoknots_2.toVienna(18, -1),
                         '(.......[{.)]}(.).')
        self.assertEqual(self.two_pseudoknots_2.toVienna(18, -1).toPairs(),
                         self.two_pseudoknots_2)
        
        self.assertEqual(self.two_pseudoknots_3.toVienna(120, -1),
                         '(.......[{.)].....}(..............................'+\
                         '..................................................'+\
                         '................)...')
        self.assertEqual(self.two_pseudoknots_3.toVienna(120, -1).toPairs(),
                         self.two_pseudoknots_3)
        
        self.assertEqual(self.three_pseudoknots.toVienna(20, -1),
                         '(....(..[)[)]([)]..]')
                         
        self.assertEqual(self.three_pseudoknots.toVienna(20, -1).toPairs(),
                         self.three_pseudoknots)
        
        self.assertEqual(self.bplist_with_one_long_pseudoknot.toVienna(45, -1),
                         '.(((((((..[[[[[[[)))))))..............]]]]]]]')
        self.assertEqual(self.bplist_with_one_long_pseudoknot.\
                         toVienna(45, -1).toPairs(),
                         self.bplist_with_one_long_pseudoknot)
        
        self.assertEqual(\
            self.bplist_with_one_gapped_pseudoknot.toVienna(45, -1),
            '.(((((((..[[[.[[[)))))))..............]]].]]]')
        self.assertEqual(\
            self.bplist_with_one_gapped_pseudoknot.toVienna(45, -1).toPairs(),
            self.bplist_with_one_gapped_pseudoknot)

        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_1.toVienna(46, -1),
            '.(((((((..[[[..))))))).]]]..(((..[[[..)))..]]]')
        
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_1.toVienna(46, -1).toPairs(),
            self.bplist_with_two_long_pseudoknots_1)
        
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_2.toVienna(46, -1),
            '.(((((((..[[[..))))))).(((..[[[..)))..]]]..]]]')
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_2.toVienna(46, -1).toPairs(),
            self.bplist_with_two_long_pseudoknots_2)
        
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_3.toVienna(46, -1),
            '.(((((((..[[[......)))))))...{{{..]]].}}}.....')
        self.assertEqual(\
            self.bplist_with_two_long_pseudoknots_3.toVienna(46, -1).toPairs(),
            self.bplist_with_two_long_pseudoknots_3)
        
        self.assertEqual(\
            self.bplist_with_two_gapped_pseudoknots_1.toVienna(46, -1),
            '.(((((((..[.[..))))))).(((..[[[..)))..]]]..].]')
        self.assertEqual(self.bplist_with_two_gapped_pseudoknots_1.\
                         toVienna(46, -1).toPairs(),
                         self.bplist_with_two_gapped_pseudoknots_1)
        
        self.assertEqual(\
            self.bplist_with_two_gapped_pseudoknots_1.toVienna(46, -1),
            '.(((((((..[.[..))))))).(((..[[[..)))..]]]..].]')
        self.assertEqual(self.bplist_with_two_gapped_pseudoknots_1.\
                         toVienna(46, -1).toPairs(),
                         self.bplist_with_two_gapped_pseudoknots_1)
        
        self.assertEqual(self.bplist_with_two_gapped_pseudoknots_2.\
                         toVienna(46, -1),
                         '.(((.(((..[[[..))).))).(.(..[[[..).)..]]]..]]]')
        self.assertEqual(self.bplist_with_two_gapped_pseudoknots_2.\
                         toVienna(46, -1).toPairs(),
                         self.bplist_with_two_gapped_pseudoknots_2)
        
    def test_toVienna_conflict(self):
        """BasePairs.toVienna should raise ConflictInBasePairsError"""
        self.assertRaises(ConflictInBasePairsError, \
            self.bplist_with_conflicts.toVienna, 100)
        self.assertRaises(ConflictInBasePairsError, \
            self.bplist_with_conflicts.toVienna, 100, -2)
        self.assertRaises(ConflictInBasePairsError, \
            BasePairs([(1,2),(2,3)]).toVienna, 4)
        
    def test_toVienna_toPairs(self):
        """BasePairs.toVienna.toPairs() should generate the same BasePairs
        """
        bps = BasePairs(((1,3), (4,5), (7,12)))
        self.assertEqual(bps.toVienna(15).toPairs(), bps)

        bps = BasePairs(((1,3), (4,5), (7,121)))
        self.assertEqual(bps.toVienna(150).toPairs(), bps)
    
    def test_toVienna_toPairs_pseudoknot(self):
        """BasePairs.toVienna.toPairs() should generate the same pseu. BasePairs 
        """
        bps = BasePairs([(1,3), (2, 5), (7,12)])
        self.assertEqual(bps.toVienna(45).toPairs(), bps)
        
    def test_solve_conflicts_of_conflicting_base_pairs(self):
        """BasePairs.make_non_conflicting_bps - single conflict
        """
        combinations = self.conflicting_bps.make_non_conflicting_bps()
        self.assertEqual(combinations, self.ref_combinations)
        
    def test_solve_conflicts_even_more_conflicting_base_pairs(self):
        """BasePairs.make_non_conflicting_bps - nested conflict
        """
        combinations = self.even_more_conflicting_bps.make_non_conflicting_bps()
        self.assertEqual(combinations, self.even_more_ref_combinations)
        
    def test_solve_conflicts_no_conflicts(self):
        """BasePairs.make_non_conflicting_bps - no conflicts in input
        """
        combinations = BasePairs([(1,2),(3,4),(5,6)]).make_non_conflicting_bps()
        self.assertEqual(combinations, [[(1,2),(3,4),(5,6)], []])
        
    def test_make_non_conflicting_viennas_of_conflicting_base_pairs(self):
        """BasePairs.make_non_conflicting_viennas - single conflict
        """
        viennas = self.conflicting_bps.make_non_conflicting_viennas(50)
        self.assertEqual(len(viennas), 5)
        self.assertEqual(viennas,
                         ('(((((.(((.((.((((.....))..)).)).))).))))).........',
                          '....................(............)................',
                          '.....................(..........).................',
                          '......................(........)..................',
                          '...................(..............)...............'))
        
    def test_make_non_conflicting_viennas_evenmore_conflicting_base_pairs(self):
        """BasePairs.make_non_conflicting_viennas - nested conflict
        """
        viennas = \
            self.even_more_conflicting_bps.make_non_conflicting_viennas(50)
        self.assertEqual(len(viennas), 6)
        self.assertEqual(viennas,
                         ('(((((.(((.((.((((.....))..)).)).))).))))).........',
                          '.....................(............)...............',
                          '....................(............)................',
                          '......................(........)..................',
                          '.....................(..........).................',
                          '...................(..............)...............'))
        
    def test_make_non_conflicting_viennas_no_conflicts(self):
        """BasePairs.make_non_conflicting_viennas - no conflicts in input
        """
        viennas = \
            BasePairs([(1,2),(3,4),(5,6)]).make_non_conflicting_viennas(10)
        self.assertEqual(len(viennas), 1)
        self.assertEqual(viennas[0], '()()()....')
        
    ###########################################################################
    ##                                                                       ##
    ##        THESE TESTS WERE TAKEN FROM PYCOGENT                           ##
    ##        (ONLY SOME MINOR CHANGES WERE MADE)                            ##
    ##                                                                       ##
    ###########################################################################
            
    def test_hasConflicts(self):
        """BasePairs.hasConflicts should return True if conflicts exist"""
        self.assertFalse(BasePairs([]).hasConflicts())
        self.assertFalse(BasePairs([(1,2),(3,4)]).hasConflicts())
        self.assertTrue(BasePairs([(1,2),(2,3)]).hasConflicts())
        self.assertTrue(BasePairs([(1,2),(2,None)]).hasConflicts())
        self.assertTrue(self.bplist_with_conflicts.hasConflicts())

    def test_directed(self):
        """BasePairs.directed should change all pairs so that a<b in (a,b)"""
        self.assertEqual(BasePairs([]).directed(),[])
        self.assertEqual(BasePairs([(2,1),(6,4),(1,7),(8,3)]).directed(),\
            BasePairs([(1,2),(1,7),(3,8),(4,6)]))
        self.assertEqual\
            (BasePairs([(5,None),(None,3)]).directed(), BasePairs([]))
        self.assertEqual(\
            BasePairs([(2,1),(1,2)]).directed(), BasePairs([(1,2)]))
            
    def test_symmetric(self):
        """BasePairs.symmetric should add (down,up) for each (up,down)"""
        self.assertEqual(BasePairs([]).symmetric(),[])
        for item in BasePairs([(1,2)]).symmetric():
            self.assertTrue(item in [(2,1),(1,2)])
        for item in BasePairs([(1,2),(1,2)]).symmetric():
            self.assertTrue(item in [(1,2),(2,1)])
        for item in BasePairs([(1,2),(3,4)]).symmetric():
            self.assertTrue(item in [(1,2),(2,1),(3,4),(4,3)])
        for item in BasePairs([(1,None)]).symmetric():
            self.assertTrue(item in [])
            
    def test_mismatches(self):
        """BasePairs.mismatches should return base pairs that can't be made"""
        # with plain string
        self.assertEqual(BasePairs([(0,1)]).mismatches('AC',{}),1)
        self.assertEqual(\
            BasePairs([(0,1)]).mismatches('AC',{('A','C'):None}),0)
        self.assertEqual(\
            BasePairs([(0,1)]).mismatches('AC',{('A','G'):None}),1)
        self.assertEqual(BasePairs([(0,1),(2,3),(3,1)]).\
        mismatches('ACGU',{('A','U'):None}),3)

        # using sequence with alphabet
        self.assertEqual(\
            BasePairs([(0,1),(0,4),(0,3)]).mismatches(Rna('ACGUA')),2)
    
    ###########################################################################
    ##                                                                       ##
    ###########################################################################

    def setUp(self):
        """Sets up environment for tests"""
        self.no_pseudoknots = BasePairs([(1,12), (2,9), (13,15)])
        self.one_pseudoknot = BasePairs([(1,12),(11,13),(14,16)])
        self.two_pseudoknots_1 = BasePairs([(1,12),(11,13),(14,16),(15,17)])
        self.two_pseudoknots_2 = BasePairs([(1,12),(9,13),(10,14),(15,17)])
        self.two_pseudoknots_3 = BasePairs([(1,12),(9,13),(10,19),(20,117)])
        self.three_pseudoknots = BasePairs(\
            [(1,12), (6,10), (9,20), (11,13), (14,16), (15,17)])
        # e.g.: (....(..{)[)]([)]..}......................
        
        self.bplist_with_one_long_pseudoknot = BasePairs([(2, 24), (3, 23), \
        (4, 22), (5, 21), (6, 20), (7, 19), (8, 18), (11, 45), (12, 44), \
        (13, 43), (14, 42), (15, 41), (16, 40), (17, 39)])
        # e.g.: .(((((((..[[[[[[[)))))))..............]]]]]]].
        
        self.bplist_with_one_gapped_pseudoknot = BasePairs([(2, 24), (3, 23), \
        (4, 22), (5, 21), (6, 20), (7, 19), (8, 18), (11, 45), (12, 44), \
        (13, 43), (15, 41), (16, 40), (17, 39)])
        # e.g.: .(((((((..[[[.[[[)))))))..............]]].]]].
        
        self.bplist_with_two_long_pseudoknots_1 = BasePairs([(2, 22), (3, 21),\
        (4, 20), (5, 19), (6, 18), (7, 17), (8, 16), (11, 26), (12, 25), \
        (13, 24), (29, 41), (30, 40), (31, 39), (34, 46), (35, 45), (36, 44)])
        # e.g.: .(((((((..[[[..))))))).]]]..(((..[[[..)))..]]]
        
        self.bplist_with_two_long_pseudoknots_2 = BasePairs([(2, 22), (3, 21),\
        (4, 20), (5, 19), (6, 18), (7, 17), (8, 16), (11, 46), (12, 45), \
        (13, 44), (24, 36), (25, 35), (26, 34), (29, 41), (30, 40), (31, 39)])
        # e.g.: .(((((((..[[[..))))))).{{{..(((..}}}..)))..]]]
        
        self.bplist_with_two_long_pseudoknots_3 = BasePairs([(2, 26), (3, 25),\
        (4, 24), (5, 23), (6, 22), (7, 21), (8, 20), (11, 37), (12, 36), \
        (13, 35), (30, 41), (31, 40), (32, 39)])
        # e.g.: .(((((((..[[[......)))))))...{{{..]]].}}}.....
        
        self.bplist_with_two_gapped_pseudoknots_1 = BasePairs([(2, 22), \
        (3, 21), (4, 20), (5, 19), (6, 18), (7, 17), (8, 16), (11, 46), \
        (13, 44), (24, 36), (25, 35), (26, 34), (29, 41), (30, 40), (31, 39)])
        # e.g.: .(((((((..[.[..))))))).{{{..(((..}}}..)))..].]
        
        self.bplist_with_two_gapped_pseudoknots_2 = BasePairs([(2, 22), \
        (3, 21), (4, 20), (6, 18), (7, 17), (8, 16), (11, 46), (12, 45), \
        (13, 44), (24, 36), (26, 34), (29, 41), (30, 40), (31, 39)])
        # e.g.: .(((.(((..[[[..))).))).(.(..[[[..).)..]]]..]]]
        
        self.bplist_with_conflicts = BasePairs([(1, 77), (1, 73), (2, 72),\
        (3, 71), (4, 70), (5, 69), (6, 68), (10, 48), (11, 47), (12, 46),\
        (17, 42), (18, 41), (19, 40), (20, 39), (21, 38), (22, 37), (23, 35),\
        (24, 34), (25, 33), (26, 32), (50, 66), (51, 65), (52, 64), (53, 63),\
        (54, 62)])
        
        self.conflicting_bps = BasePairs([(1, 41), (2, 40), (3, 39), (4, 38),
            (5, 37), (7, 35), (8, 34), (9, 33), (11, 31), (12, 30), (14, 28),
            (15, 27), (16, 24), (17, 23), (20, 35), (21, 34), (22, 33),
            (23, 32)])
        
        # two lists should be created
        self.ref_combinations = [[(1, 41), (2, 40), (3, 39), (4, 38), (5, 37),
            (7, 35), (8, 34), (9, 33), (11, 31), (12, 30), (14, 28), (15, 27),
            (16, 24), (17, 23)], [(21, 34), (22, 33), (23, 32), (20, 35)]]
        
        self.even_more_conflicting_bps = \
            BasePairs(self.conflicting_bps + [(22, 35)])
        self.even_more_ref_combinations = [
            [(1, 41), (2, 40), (3, 39), (4, 38), (5, 37), (7, 35), (8, 34),
                (9, 33), (11, 31), (12, 30), (14, 28), (15, 27), (16, 24),
                (17, 23)],
            [(22, 35), (21, 34), (23, 32)],
            [(22, 33), (20, 35)]
            ]


class SolveConflictsTests(TestCase):
    """tests for solve_conflicts 
    """
    def setUp(self):
        self.conflicting_bps = BasePairs([(1, 41), (2, 40), (3, 39), (4, 38),
            (5, 37), (7, 35), (8, 34), (9, 33), (11, 31), (12, 30), (14, 28),
            (15, 27), (16, 24), (17, 23), (20, 35), (21, 34), (22, 33),
            (23, 32)])
        
        # two lists should be created
        self.ref_combinations = [[(1, 41), (2, 40), (3, 39), (4, 38), (5, 37),
            (7, 35), (8, 34), (9, 33), (11, 31), (12, 30), (14, 28), (15, 27),
            (16, 24), (17, 23)], [(21, 34), (22, 33), (23, 32), (20, 35)]]
        
        self.even_more_conflicting_bps = \
            BasePairs(self.conflicting_bps + [(22, 35)])
        self.even_more_ref_combinations = [
            [(1, 41), (2, 40), (3, 39), (4, 38), (5, 37), (7, 35), (8, 34),
                (9, 33), (11, 31), (12, 30), (14, 28), (15, 27), (16, 24),
                (17, 23)],
            [(22, 35), (21, 34), (23, 32)],
            [(22, 33), (20, 35)]
            ]
        
    def test_solve_conflicts_of_conflicting_base_pairs(self):
        """solve_conflicts - single conflict
        """
        combinations = solve_conflicts(self.conflicting_bps)
        self.assertEqual(combinations, self.ref_combinations)
        
        
    def test_solve_conflicts_even_more_conflicting_base_pairs(self):
        """solve_conflicts - nested conflict
        """
        combinations = solve_conflicts(self.even_more_conflicting_bps)
        self.assertEqual(combinations, self.even_more_ref_combinations)
        
    def test_solve_conflicts_no_conflicts(self):
        """solve_conflicts - no conflicts in input
        """
        combinations = solve_conflicts(BasePairs([(1,2),(3,4),(5,6)]))
        self.assertEqual(combinations, [[(1,2),(3,4),(5,6)], []])

if  __name__ == "__main__":
    main()

