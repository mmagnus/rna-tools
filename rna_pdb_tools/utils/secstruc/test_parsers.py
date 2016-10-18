#!/usr/bin/env python
#-*-coding: utf-8-*-

"""
tests for parsers.py.
"""

__author__ = "Tomasz Puton"
__credits__ = "Ewa Tkalińska, Łukasz Kozłowski, Kristian Rother"
__license__ = "GNU GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

import os, tempfile

TEST_DIR_PATH = os.path.join(os.path.abspath('test_data'), 'CRW_00583.ct')

from unittest import TestCase,  main

from parsers import RNAStrandRecord, parse_hotknots, \
parse_afold, has_overlapping_basepairs, remove_pseudoknot, is_bp_a_pseudoknot,\
parse_ct, parse_pknots, parse_multiple_vienna, parse_sfold, parse_bp, \
parse_vienna, parse_rnashapes, make_dotbracket_from_bplist, \
make_bplist_from_dotbracket, PseudoknotError, compose_ss_from_cmfinder_motives
from secstruc import PseudoknotTokensError, BasePairs, ViennaStructure
 
class RNASecstrucParserTests(TestCase):
    
    def test_parse_afold(self):
        """Should create a bplist from Afold output data (line iterator)"""
        afold_out = '''..G..1..77..A.. = -28.20 Multidomain 
..G..1..73..C.. = -26.50 0x0 
..G..2..72..C.. = -23.20 0x0 
..A..3..71..U.. = -20.80 0x0 
..G..4..70..U.. = -20.20 0x0 
..C..5..69..G.. = -17.70 0x0 
..G..6..68..C.. = -15.30 0x0 
..G..7..67..C.. = -12.00 multi 
..A..27..48..U.. = -9.40 0x0 
..C..28..47..G.. = -7.20 0x0 
..C..29..46..G.. = -3.90 0x0 
..U..30..45..G.. = -1.80 1x0 
..C..32..44..G.. = -4.10 0x0 
..C..33..43..G.. = -0.80 0x0 
..U..34..42..A.. = 1.30 0x0 
..G..35..41..C.. = 3.40 hairpin 
..G..50..66..U.. = -7.10 0x0 
..C..51..65..G.. = -4.60 0x0 
..G..52..64..C.. = -2.20 0x0 
..G..53..63..C.. = 1.10 0x0 
..G..54..62..C.. = 4.40 hairpin
'''.split('\n')
        result = parse_afold(afold_out)
        self.assertEqual(result[0], [(1, 77), (1, 73), (2, 72), (3, 71),
            (4, 70), (5, 69), (6, 68), (7, 67), (27, 48), (28, 47), (29, 46),
            (30, 45), (32, 44), (33, 43), (34, 42), (35, 41), (50, 66),
            (51, 65), (52, 64), (53, 63), (54, 62)])
        self.assertEqual(result[1], -28.20)
        
    def test_sfold(self):
        """Should create a bplist from Sfold output data (line iterator)."""
        sfold_out = """Structure 1
1 20
2 19
3 18
4 17
5 16
6 15
7 14
""".split('\n')
        result = parse_sfold(sfold_out)
        self.assertEqual(result, [(1, 20), (2, 19), (3, 18), (4, 17), (5, 16), (6, 15), (7, 14)])


    def test_parse_pknots(self):
        """Should parse a bplist from pknots output data (line iterator)"""
        pknots_data = iter("""NAM  fasta
SEQ  +SS
       1 G   A   U   C   A   U   A   U   G   A   C   U   G   A   A   G   A   G   G   C
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  
         .   .   .   .   .   .   .   .   .   .   18  17  .   .   .   .   .   11  10  .    

      21 A   U
         20  21  
         21  20  
""".split('\n'))
        result = parse_pknots(pknots_data)
        self.assertEqual(result[0], [(11, 19), (12, 18), (21, 22)])
        self.assertEqual(result[1], None)
        
    def test_parse_pknots_tricky(self):
        """ """
        pknots_data = iter("""NAM  input
DES  seq
SEQ  +SS
       1 A   A   A   A   A   A   A   A   A   A   A   A   A   A   A   A   A   A   A   A
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  
         20  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   

      21 A
         20  
         0
         
----------------------------------------------------------------------
   Log odds score:    -9999999 
energy (kcal/mol):    10000.00 
""".split('\n'))
        result = parse_pknots(pknots_data)
        self.assertEqual(result[0], BasePairs([(1, 21)]))
        self.assertAlmostEqual(result[1], 10000.00, places=2)

    def test_parse_ct(self):
        """Should parse a bplist from ct data (line iterator)"""
        ct_data = """73 ENERGY =     -17.50    S.cerevisiae_tRNA-PHE
 1 G       0    2   72    1
 2 C       1    3   71    2
 3 G       2    4   70    3
 """.split('\n')
        result = parse_ct(ct_data)
        self.assertEqual(result, [(1, 72), (2, 71), (3, 70)])
        
        
    def test_parse_bp(self):
        """Should parse bp files with and without header from a line iterator."""
        bp_with_header = """Filename: test.bpseq
Organism: Some organism
Accession Number: XYZ123
Citation and related information available at http://www.rna.ccbb.utexax.edu
1 C 0
2 C 9
3 U 8
4 G 0
5 A 0
6 A 0
7 C 0
8 A 3
9 G 2
""".split('\n')
        result = parse_bp(bp_with_header)
        self.assertEqual(result, [(2, 9), (3, 8)])

        bp_without_header = """1 C 0
2 C 9
3 U 8
4 G 0
5 A 0
6 A 0
7 C 0
8 A 3
9 G 2
""".split('\n')
        result_without = parse_bp(bp_without_header)
        self.assertEqual(result_without, result)

    def test_parse_bp_conflict(self):
        """parse_bp should work for conflicting base pairs
        """
        input_bpseq =  """1 C 0
2 C 9
3 U 8
3 U 9
4 G 0
5 A 0
6 A 0
7 C 0
8 A 3
9 G 2
9 G 3
"""
        result = parse_bp(input_bpseq)
        self.assertEqual(result, [(2, 9), (3, 8), (3,9)])
        
        fd, path = tempfile.mkstemp()
        open(path, 'w').write(input_bpseq)
        
        result = parse_bp(open(path))
        self.assertEqual(result, [(2, 9), (3, 8), (3,9)])
        os.close(fd)
        os.unlink(path)
        

    def test_parse_vienna_no_energy(self):
        """Should parse a secstruc from vienna data without energy"""
        vienna_data = """> blabla name
AAAGUUCCCAAAUUU
(((.((...)).)))""".split('\n')
        result = parse_vienna(vienna_data)
        self.assertEqual(result, \
            (' blabla name', 'AAAGUUCCCAAAUUU', '(((.((...)).)))', None))

    def test_parse_vienna_with_energy(self):
        """Should parse a secstruc from vienna data with energy"""
        vienna_data = """> blabla name
AAAGUUCCCAAAUUU
(((.((...)).))) (-12.123)""".split('\n')
        result = parse_vienna(vienna_data)
        self.assertEqual(result, \
            (' blabla name', 'AAAGUUCCCAAAUUU', '(((.((...)).)))', -12.123))
        
    def test_parse_rnashapes(self):
        """Should correctly parse RNAshapes data."""
        rnashapes_data = """>fasta
        GATCATATGACTGAAGAGGCATATATGACTGAAGAGGCATATATGACAAGATGACTGAA
-8.50   .(((......((....)).((((((((.((.....))))))))))....))).......  [[][]]
-8.30   ...................((((((((.((.....))))))))))..............  []
""".split('\n')
        result = parse_rnashapes(rnashapes_data)
        self.assertEqual(result.next(), \
            (-8.5, '.(((......((....)).((((((((.((.....))))))))))....))).......', '[[][]]'))
        self.assertEqual(result.next(), \
            (-8.3, '...................((((((((.((.....))))))))))..............', '[]'))
        
        rnashapes_data2 = """        GATCATATGACTGAAGAGGCATATATGACTGAAGAGGCATATATGACAAGATGACTGAA
-8.50   .(((......((....)).((((((((.((.....))))))))))....))).......  [[][]]
-8.40   .(((......((....)).((((((((.((....)).))))))))....))).......  [[][]]
-8.30   ...................((((((((.((.....))))))))))..............  []
-8.20   ...................((((((((.((....)).))))))))..............  []
-9.00   ..((((....((....)).((((((((.((.....)))))))))).....)))).....  [[][]]
""".split('\n')
        result = [r for r in parse_rnashapes(rnashapes_data2)]
        self.assertEqual(len(result), 5)
        
    def test_parse_multiple_vienna(self):
        """Should yield secondary structure strings from a multiple vienna line iterator."""
        multiple_vienna_data = """> fasta [100]
UAUGACUAUA   -920    100
.((....)).  -8.40
(((....)))  -8.50
""".split('\n')
        result = parse_multiple_vienna(multiple_vienna_data)
        self.assertEqual(result.next(), (None, None, '.((....)).', -8.4))
        self.assertEqual(result.next(), (None, None, '(((....)))', -8.5))
        
##    def test_get_nr_of_pseudoknots(self):
##        """Should return number of pseudoknots in a given bplist"""
##        bplist_with_one_pseudoknot = [[2, 10], [4, 6], [3, 9], [5, 8]]
##        bplist_with_two_pseudoknots = [[1, 10], [3, 9], [5, 7], [2, 8], [4, 6]]
##        bplist_without_any_pseudoknots = [[5, 8], [2, 10], [3, 9]]
##        bplist1 = [[1, 4], [2, 5], [6, 9]]
##        bplist0 = [[2, 5], [6, 9]]
##        bplist3 = [[1,4], [2,5], [6,9], [120, 140], [130, 144], [18, 20], \
##        [19, 21]]
##        bplist4 = [[1,4], [2,5], [6,9], [120, 140], [121, 170], [130, 144], \
##        [18, 20], [19, 21]]
##        
##        bplist_with_one_long_pseudoknot = [[2, 24], [3, 23], [4, 22], \
##        [5, 21], [6, 20], [7, 19], [8, 18], [11, 45], [12, 44], [13, 43], \
##        [14, 42],\
##        [15, 41], [16, 40], [17, 39]]
##        # e.g.: .(((((((..[[[[[[[)))))))..............]]]]]]].
##        
##        bplist_with_one_gapped_pseudoknot = [[2, 24], [3, 23], [4, 22],\
##        [5, 21], [6, 20], [7, 19], [8, 18], [11, 45], [12, 44], [13, 43],\
##        [15, 41], [16, 40], [17, 39]]
##        # e.g.: .(((((((..[[[.[[[)))))))..............]]].]]].
##        
##        bplist_with_two_long_pseudoknots_1 = [[2, 22], [3, 21], [4, 20],\
##        [5, 19], [6, 18], [7, 17], [8, 16], [11, 26], [12, 25], [13, 24],\
##        [29, 41], [30, 40], [31, 39], [34, 46], [35, 45], [36, 44]]
##        # e.g.: .(((((((..[[[..))))))).]]]..(((..[[[..)))..]]]
##        
##        bplist_with_two_long_pseudoknots_2 = [[2, 22], [3, 21], [4, 20],\
##        [5, 19], [6, 18], [7, 17], [8, 16], [11, 26], [12, 25], [13, 24],\
##        [29, 41], [30, 40], [31, 39], [34, 46], [35, 45], [36, 44]]
##        # e.g.: .(((((((..[[[..))))))).]]]..(((..{{{..)))..}}}
##        
##        bplist_with_two_long_pseudoknots_3 = [[2, 22], [3, 21], [4, 20],\
##        [5, 19], [6, 18], [7, 17], [8, 16], [11, 46], [12, 45], [13, 44],\
##        [24, 36], [25, 35], [26, 34], [29, 41], [30, 40], [31, 39]]
##        # e.g.: .(((((((..[[[..))))))).{{{..(((..}}}..)))..]]]
##        
##        bplist_with_two_long_pseudoknots_4 = [[2, 22], [3, 21], [4, 20],\
##        [5, 19], [6, 18], [7, 17], [8, 16], [11, 46], [12, 45], [13, 44],\
##        [24, 36], [25, 35], [26, 34], [29, 41], [30, 40], [31, 39]]
##        # e.g.: .(((((((..xxx..))))))).{{{..(((..}}}..)))..XXX
##        
##        bplist_with_two_long_pseudoknots_5 = [[2, 26], [3, 25], [4, 24],\
##        [5, 23], [6, 22], [7, 21], [8, 20], [11, 37], [12, 36], [13, 35],\
##        [30, 41], [31, 40], [32, 39]]
##        # e.g.: .(((((((..[[[......)))))))...{{{..]]].}}}.....
##        
##        bplist_with_two_long_pseudoknots_6 = [[2, 26], [3, 25], [4, 24],\
##        [5, 23], [6, 22], [7, 21], [8, 20], [11, 37], [12, 36], [13, 35],\
##        [30, 41], [31, 40], [32, 39]]
##        # e.g.: .(((((((..[[[......)))))))...(((..]]].))).....
##        
##        bplist_with_two_gapped_pseudoknots_1 = [[2, 22], [3, 21], [4, 20],\
##        [5, 19], [6, 18], [7, 17], [8, 16], [11, 46], [13, 44], [24, 36],\
##        [25, 35], [26, 34], [29, 41], [30, 40], [31, 39]]
##        # e.g.: .(((((((..[.[..))))))).{{{..(((..}}}..)))..].]
##        
##        bplist_with_two_gapped_pseudoknots_2 = [[2, 22], [3, 21], [4, 20],\
##        [6, 18], [7, 17], [8, 16], [11, 46], [12, 45], [13, 44], [24, 36],\
##        [26, 34], [29, 41], [30, 40], [31, 39]]
##        # e.g.: .(((.(((..[[[..))).))).{.{..(((..}.}..)))..]]]
##        
##        bplist_with_three_pseudoknots = [[2, 42], [3, 41], [4, 40], [5, 39],\
##        [6, 38], [7, 37], [8, 36], [10, 31], [11, 27], [12, 26], [13, 25],\
##        [18, 45], [19, 44], [20, 43]]
##        # e.g.: .(((((((.{(((....[[[....)))...}....)))))))]]].
##        
##        self.assertEqual(get_nr_of_pseudoknots(bplist_with_one_pseudoknot), 1)
##        self.assertEqual(get_nr_of_pseudoknots(bplist1), 1)
##        self.assertEqual(get_nr_of_pseudoknots(bplist_with_two_pseudoknots), 2)
##        self.assertEqual(get_nr_of_pseudoknots(bplist3), 3)
##        self.assertEqual(get_nr_of_pseudoknots(bplist4), 4)
##        self.assertEqual(\
##        get_nr_of_pseudoknots(bplist_without_any_pseudoknots), 0)
##        self.assertEqual(get_nr_of_pseudoknots(bplist0), 0)
##
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_one_long_pseudoknot), 1)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_one_gapped_pseudoknot), 1)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_long_pseudoknots_1), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_long_pseudoknots_2), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_long_pseudoknots_3), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_long_pseudoknots_4), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_long_pseudoknots_5), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_long_pseudoknots_6), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_gapped_pseudoknots_1), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_two_gapped_pseudoknots_2), 2)
##        self.assertEqual(\
##            get_nr_of_pseudoknots(bplist_with_three_pseudoknots), 3)
        
    def test_is_bp_a_pseudoknot(self):
        """Should return True if bp is pseudoknot in a given bplist, False otherwise"""
        bplist3 = [[1,4], [2,5], [6,9], [120, 140], [130, 144], [18, 20], \
        [19, 21]]
        self.assertEqual(is_bp_a_pseudoknot(1, 4, bplist3), True)
        self.assertEqual(is_bp_a_pseudoknot(120, 140, bplist3), True)
        self.assertEqual(is_bp_a_pseudoknot(6, 9, bplist3), False)
        self.assertEqual(is_bp_a_pseudoknot(19, 21, bplist3), True)
    
    def test_make_dotbracket_from_bplist(self):
        """Takes a length and list of base pairs, and returns secstruc."""
        bplist = [[2, 10], [3, 9], [5, 8]]
        result = make_dotbracket_from_bplist(10, bplist)
        self.assertEqual(result, '.((.(..)))')
        
    def test_make_dotbracket_from_bplist_pseudoknot(self):
        """Should correctly produce sec struct with pseudoknots in Vienna."""
        # one pseudonot
        bplist_1 = [[1, 4], [2, 5], [6, 9]]
        pseudoknot_tokens_1 = [('<','>')]
        result_1 = \
        make_dotbracket_from_bplist(10, bplist_1, pseudoknot_tokens_1)
        self.assertEqual(result_1, '(<.)>(..).')
        
        # three pseudoknots
        bplist_3 = [[1, 4], [2, 5], [6, 9], [7,10], [11, 13], [12, 14]]
        # default pseudoknots_tokens are [('[',']'), ('{','}'), ('<','>')]
        result_3 = make_dotbracket_from_bplist(20, bplist_3)
        self.assertEqual(result_3, '([.)]([.)]([)]......')
        
    def test_make_dotbracket_from_bplist_difficult(self):
        """Should raise PseudoknotTokensError if there's not enough tokens."""
        # three pseudoknots...
        bplist_3 = [(1, 4), (2, 5), (6, 9), (7,10), (11, 13), (12, 14)]
        # ...but just one pair of tokens
        pseudoknot_tokens_1 = [('<','>')]
        
        #self.assertRaises(PseudoknotTokensError, make_dotbracket_from_bplist,
        #                  20, bplist_3, pseudoknot_tokens_1)
        
        dbn = make_dotbracket_from_bplist(20, bplist_3, pseudoknot_tokens_1)
        self.assertEqual(dbn, '(<.)>(<.)>(<)>......')
        self.assertEqual(ViennaStructure(dbn).toPairs(), bplist_3)
        
    def test_make_bplist_from_dotbracket(self):
        """Should correctly make bplist from a dot-bracket (Vienna) notation"""
        result = make_bplist_from_dotbracket('.((.(..)))')
        self.assertEqual(result, [(2, 10), (3, 9), (5, 8)])

        result_2 = make_bplist_from_dotbracket\
                   ('.(((((((.{(((....[[[....)))...}....)))))))]]].')
        self.assertEqual(result_2, [(2, 42), (3, 41), (4, 40), (5, 39),\
        (6, 38), (7, 37), (8, 36), (10, 31), (11, 27), (12, 26), (13, 25),\
        (18, 45), (19, 44), (20, 43)])

    def test_make_bplist_from_dotbracket_pknot(self):
        """Should correctly make bplist from a Vienna containing pseudoknot"""
        result = make_bplist_from_dotbracket('.({.(..}))')
        self.assertEqual(result, [(2, 10), (3, 8), (5, 9)])
        result = make_bplist_from_dotbracket('.(c.(..C))')
        self.assertEqual(result, [(2, 10), (3, 8), (5, 9)])
        
    def test_remove_pseudoknot(self):
        """Should correctly remove pseudoknot from bplist."""
        bplist = [[2, 11], [3, 8], [4, 10], [5, 9]]
        result = remove_pseudoknot(bplist)
        self.assertEqual(result, [(2, 11), (4, 10), (5, 9)])
        
    def test_remove_pseudoknot_error(self):
        """Should raise PseudoknotError during pseudoknot removal."""
        self.assertRaises(PseudoknotError, remove_pseudoknot, [[2, 11], [3, 8], [3, 10], [5, 9]])
        
    def test_has_overlapping_basepairs(self):
        """Checks if some basepairs have identical indices"""
        self.assertEqual(has_overlapping_basepairs([[2, 11], [3, 8], [4, 10], [5, 9]]), False)
        self.assertEqual(has_overlapping_basepairs([[2, 11], [3, 8], [3, 10], [5, 9]]), True)

#    def test_separate_pseudoknot(self):
#        pseudoknot = '(((..(((.aa..))b.AA.)B..))))))'
#        bplist = make_bplist_from_dotbracket(pseudoknot)
#        knotless = ['(((..(((.....)).....)...))))))', 
#                        '.........((......))...........', 
#                        '...............(.....)........']
#        result = [k for k in separate_pseudoknot(len(pseudoknot), bplist)]
#        self.assertEqual(result,  knotless)

class HotknotsParserTests(TestCase):
    """A test suite for parse_hotknots function.
    """
    def setUp(self):
        """Sets up environment for tests.
        """
        self.data = \
                ['0 0 child nodes have been created\n',
                'In total, 112 nodes created. \n',
                ' total number of Rna structures: 11 \n',
                ' sub optimal energy  -11.850000 kal/mol\n',
                '3-5; 48-50,  6-7; 45-46,  10-16; 36-42,  22-24; 33-35,  \n',
                ' sub optimal energy  -10.650000 kal/mol\n',
                '3-7; 33-37,  9-11; 27-29,  13-15; 24-26,  39-43; 47-51,  \n',
                ' sub optimal energy  -10.480000 kal/mol\n',
                '3-7; 33-37,  8-10; 22-24,  39-42; 48-51,  \n',
                ' sub optimal energy  -9.640000 kal/mol\n',
                '4-5; 27-28,  8-10; 22-24,  35-39; 45-49,  \n',
                ' sub optimal energy  -9.200000 kal/mol\n',
                '3-7; 33-37,  9-11; 27-29,  39-43; 47-51,  \n',
                ' sub optimal energy  -8.400000 kal/mol\n',
                '8-12; 32-36,  39-42; 48-51,  \n',
                ' sub optimal energy  -8.290000 kal/mol\n',
                '3-7; 33-37,  11-15; 28-32,  39-43; 47-51,  \n',
                ' sub optimal energy  -8.200000 kal/mol\n',
                '3-5; 10-12,  22-24; 33-35,  39-43; 47-51,  \n',
                ' sub optimal energy  -8.100000 kal/mol\n',
                '10-16; 36-42,  22-24; 33-35,  \n',
                ' sub optimal energy  -7.990000 kal/mol\n',
                '7-10; 22-25,  27-29; 34-36,  39-42; 48-51,  \n',
                ' sub optimal energy  -7.890000 kal/mol\n',
                '7-10; 22-25,  35-39; 45-49,  \n', '\n']
        
        self.fd, self.path = tempfile.mkstemp()
        fh = open(self.path, 'w')
        for line in self.data:
            fh.write(line)
        fh.close()
        
        self.pseudoknot_easy = ['0 0 child nodes have been created\n', 
                           'In total, 6 nodes created. \n',
                           ' total number of Rna structures: 2 \n',
                           ' sub optimal energy  -8.300000 kal/mol\n',
                           '1-5; 12-16,  \n',
                           ' sub optimal energy  -6.722000 kal/mol\n',
                           '1-5; 12-16,  8-11; 24-27, \n']
        
        self.fd_p_easy, self.path_p_easy = tempfile.mkstemp()
        fh = open(self.path_p_easy, 'w')
        for line in self.pseudoknot_easy:
            fh.write(line)
        fh.close()
        
        self.pseudoknot_medium = ['0 0 child nodes have been created\n', 
                           'In total, 3 nodes created. \n',
                           ' total number of Rna structures: 1 \n',
                           ' sub optimal energy  -8.300000 kal/mol\n',
                           '1-5; 12-16, 8-9; 27-28, 19-25; 29-35, \n'
                           ' sub optimal energy  999999 kal/mol\n',
                           '1-5; 12-16, 8-9; 27-28, 19-25; 29-35, 36-37; 39-40, \n']
        
        self.fd_p_medium, self.path_p_medium = tempfile.mkstemp()
        fh = open(self.path_p_medium, 'w')
        for line in self.pseudoknot_medium:
            fh.write(line)
        fh.close()

        self.pseudoknot_hard = ['0 0 child nodes have been created\n', 
                           'In total, 3 nodes created. \n',
                           ' total number of Rna structures: 1 \n',
                           ' sub optimal energy  -7.123456 kal/mol\n',
                           '1-5; 12-16, 8-9; 27-28, 10-11; 39-40, 19-25; 29-35, \n']
        
        self.fd_p_hard, self.path_p_hard = tempfile.mkstemp()
        fh = open(self.path_p_hard, 'w')
        for line in self.pseudoknot_hard:
            fh.write(line)
        fh.close()
        
    def tearDown(self):
        """Removes leftovers.
        """
        os.close(self.fd)
        os.unlink(self.path)

        os.close(self.fd_p_easy)
        os.unlink(self.path_p_easy)

        os.close(self.fd_p_medium)
        os.unlink(self.path_p_medium)

        os.close(self.fd_p_hard)
        os.unlink(self.path_p_hard)
        
    def test_parser(self):
        """parse_hotknots should work correctly on structs without pseudkonts"""
        parser = parse_hotknots(open(self.path), 53)
        res1 = parser.next()
        self.assertEqual(res1[0], \
                        '..(((((..(((((((.....(((........))))))))))..)).)))...')
        self.assertEqual(res1[1], -11.850000)
        
        res2 = parser.next()
        self.assertEqual(res2[0], \
                        '..(((((.(((.(((........))))))...))))).(((((...)))))..')
        self.assertEqual(res2[1], -10.650000)

        res3 = parser.next()
        self.assertEqual(res3[0], \
                        '..((((((((...........)))........))))).((((.....))))..')
        self.assertEqual(res3[1], -10.480000)

        res4 = parser.next()
        self.assertEqual(res4[0], \
                        '...((..(((...........)))..))......(((((.....)))))....')
        self.assertEqual(res4[1], -9.640000)
        
        res5 = parser.next()
        self.assertEqual(res5[0], \
                        '..(((((.(((...............)))...))))).(((((...)))))..')
        self.assertEqual(res5[1], -9.200000)

        res6 = parser.next()
        self.assertEqual(res6[0], \
                        '.......(((((...................)))))..((((.....))))..')
        self.assertEqual(res6[1], -8.400000)

        res7 = parser.next()
        self.assertEqual(res7[0], \
                        '..(((((...(((((............)))))))))).(((((...)))))..')
        self.assertEqual(res7[1], -8.290000)

        res8 = parser.next()
        self.assertEqual(res8[0], \
                        '..(((....))).........(((........)))...(((((...)))))..')
        self.assertEqual(res8[1], -8.200000)
        
        res9 = parser.next()
        self.assertEqual(res9[0], \
                        '.........(((((((.....(((........))))))))))...........')
        self.assertEqual(res9[1], -8.100000)

        res10 = parser.next()
        self.assertEqual(res10[0], \
                        '......((((...........)))).(((....)))..((((.....))))..')
        self.assertEqual(res10[1], -7.990000)

        res11 = parser.next()
        self.assertEqual(res11[0], \
                        '......((((...........)))).........(((((.....)))))....')
        self.assertEqual(res11[1], -7.890000)
        
        self.assertRaises(StopIteration, parser.next)
        
    def test_parser_pseudoknot_easy(self):
        """parse_hotknots should work correctly on pseudoknotted structs, p.1"""
        parser = parse_hotknots(open(self.path_p_easy), 27)
        res1 = parser.next()
        self.assertEqual(res1[0], '(((((......)))))...........')
        self.assertEqual(res1[1], -8.300000)
        
        res2 = parser.next()
        self.assertEqual(res2[0], '(((((..[[[[))))).......]]]]')
        self.assertEqual(res2[1], -6.722000)
        
        self.assertRaises(StopIteration, parser.next)
        
    def test_parser_pseudoknot_medium(self):
        """parse_hotknots should work correctly on pseudoknotted structs, p.2"""
        parser = parse_hotknots(open(self.path_p_medium), 40)
        res1 = parser.next()
        self.assertEqual(res1[0], '(((((..[[..)))))..(((((((.]]))))))).....')
        self.assertEqual(res1[1], -8.300000)

        res2 = parser.next()
        self.assertEqual(res2[0], '(((((..[[..)))))..(((((((.]])))))))((.))')
        self.assertEqual(res2[1], 999999)
        
        self.assertRaises(StopIteration, parser.next)
        
    def test_parser_pseudoknot_hard(self):
        """parse_hotknots should work correctly on pseudoknotted structs, p.3"""
        parser = parse_hotknots(open(self.path_p_hard), 40)
        res1 = parser.next()
        self.assertEqual(res1[0], '(((((..[[{{)))))..(((((((.]])))))))...}}')
        self.assertEqual(res1[1], -7.123456)
        
        self.assertRaises(StopIteration, parser.next)


class RNAStrandRecordTests(TestCase):

    def test_parse(self):
        """RNAStrandRecord._parse should work correctly"""
        defline, sequence, structure, pairs = self.record._parse()
        
        self.assertEqual(defline, self.ref_defline)
        self.assertEqual(sequence, self.ref_sequence)
        self.assertEqual(structure, self.ref_structure)
        
        self.assertEqual(pairs, self.ref_pairs)
        
    def test_vienna(self):
        """RNAStrandRecord.vienna proporty should work correctly"""
        self.assertEqual(self.record.vienna, self.ref_vienna)
        
    def test_is_valid_sequence(self):
        """RNAStrandRecord._is_valid_sequence should work correctly"""
        self.assertEqual(self.record._is_valid_sequence('ACUAGCUACGCCCC'), True)
        self.assertEqual(self.record._is_valid_sequence('ACUAGCUACGCNC'), False)
        self.assertEqual(self.record._is_valid_sequence('ACU#GCUA*GC7C'), True)
        self.assertEqual(self.record._is_valid_sequence('ACU#GCUA~GC7C'), False)
        
    def test_rec_not_enough_tokens(self):
        """RNAStrandRecord should deal with highly pseudoknotted RNA structure
        """
        record = RNAStrandRecord(self.not_enough_pseudoknot_tokens_path)
        
        record_pairs = record.pairs
        structure_pairs = record.structure.toPairs()
        
        self.assertEqual(record.structure.toPairs(), record.pairs)
        self.assertEqual(record.pairs, self.ref_not_en_tokens_pairs)
        self.assertEqual(record.valid, True)
        self.assertEqual(len(record.sequence), 2930)
        
    def setUp(self):
        """ """
        self.record_str = """# File ASE_00001.ct
# RNA SSTRAND database, RNase P RNA, Acidianus ambivalens, strain DSM 3772
# External source: RNase P Database, file name: A.ambivalens.ct2
  262 ENERGY = 0 A.ambivalens RNase P RNA
    1 g       0    2    0    1
    2 a       1    3    0    2
    3 g       2    4    0    3
    4 g       3    5    0    4
    5 a       4    6    0    5
    6 a       5    7  261    6
    7 a       6    8  260    7
    8 g       7    9  259    8
    9 u       8   10    0    9
   10 c       9   11  258   10
   11 c      10   12  257   11
   12 c      11   13  256   12
   13 g      12   14  255   13
   14 c      13   15  254   14
   15 c      14   16  182   15
   16 U      15   17  181   16
   17 C      16   18  180   17
   18 C      17   19  179   18
   19 A      18   20    0   19
   20 G      19   21  209   20
   21 A      20   22  208   21
   22 U      21   23  207   22
   23 C      22   24  206   23
   24 A      23   25  205   24
   25 A      24   26    0   25
   26 G      25   27  178   26
   27 G      26   28  177   27
   28 G      27   29  176   28
   29 A      28   30  175   29
   30 A      29   31  174   30
   31 G      30   32   44   31
   32 U      31   33   43   32
   33 C      32   34   42   33
   34 C      33   35   41   34
   35 C      34   36   40   35
   36 G      35   37    0   36
   37 C      36   38    0   37
   38 G      37   39    0   38
   39 A      38   40    0   39
   40 G      39   41   35   40
   41 G      40   42   34   41
   42 G      41   43   33   42
   43 A      42   44   32   43
   44 C      43   45   31   44
   45 A      44   46   59   45
   46 A      45   47   58   46
   47 G      46   48   57   47
   48 G      47   49   56   48
   49 G      48   50   55   49
   50 U      49   51    0   50
   51 A      50   52    0   51
   52 G      51   53    0   52
   53 U      52   54    0   53
   54 A      53   55    0   54
   55 C      54   56   49   55
   56 C      55   57   48   56
   57 C      56   58   47   57
   58 U      57   59   46   58
   59 U      58   60   45   59
   60 G      59   61  173   60
   61 G      60   62  172   61
   62 C      61   63    0   62
   63 A      62   64    0   63
   64 A      63   65    0   64
   65 C      64   66  170   65
   66 U      65   67  169   66
   67 G      66   68  168   67
   68 C      67   69  167   68
   69 A      68   70    0   69
   70 C      69   71    0   70
   71 A      70   72    0   71
   72 G      71   73    0   72
   73 A      72   74    0   73
   74 A      73   75    0   74
   75 A      74   76    0   75
   76 A      75   77    0   76
   77 C      76   78    0   77
   78 U      77   79    0   78
   79 U      78   80    0   79
   80 A      79   81    0   80
   81 C      80   82  159   81
   82 C      81   83  158   82
   83 C      82   84  157   83
   84 C      83   85  156   84
   85 U      84   86  155   85
   86 A      85   87  154   86
   87 A      86   88  153   87
   88 A      87   89  152   88
   89 U      88   90  151   89
   90 A      89   91  150   90
   91 U      90   92  149   91
   92 U      91   93    0   92
   93 C      92   94    0   93
   94 A      93   95    0   94
   95 A      94   96    0   95
   96 U      95   97    0   96
   97 G      96   98    0   97
   98 A      97   99    0   98
   99 G      98  100    0   99
  100 G      99  101    0  100
  101 A     100  102    0  101
  102 U     101  103    0  102
  103 U     102  104    0  103
  104 U     103  105    0  104
  105 G     104  106    0  105
  106 A     105  107    0  106
  107 U     106  108    0  107
  108 U     107  109  140  108
  109 C     108  110  139  109
  110 G     109  111    0  110
  111 A     110  112    0  111
  112 C     111  113    0  112
  113 U     112  114  136  113
  114 C     113  115  135  114
  115 U     114  116  134  115
  116 U     115  117  133  116
  117 A     116  118  132  117
  118 C     117  119  131  118
  119 C     118  120  130  119
  120 U     119  121  129  120
  121 U     120  122  128  121
  122 G     121  123  127  122
  123 G     122  124    0  123
  124 C     123  125    0  124
  125 G     124  126    0  125
  126 A     125  127    0  126
  127 C     126  128  122  127
  128 A     127  129  121  128
  129 A     128  130  120  129
  130 G     129  131  119  130
  131 G     130  132  118  131
  132 U     131  133  117  132
  133 A     132  134  116  133
  134 A     133  135  115  134
  135 G     134  136  114  135
  136 A     135  137  113  136
  137 U     136  138    0  137
  138 A     137  139    0  138
  139 G     138  140  109  139
  140 A     139  141  108  140
  141 U     140  142    0  141
  142 G     141  143    0  142
  143 A     142  144    0  143
  144 A     143  145    0  144
  145 G     144  146    0  145
  146 A     145  147    0  146
  147 G     146  148    0  147
  148 A     147  149    0  148
  149 A     148  150   91  149
  150 U     149  151   90  150
  151 A     150  152   89  151
  152 U     151  153   88  152
  153 U     152  154   87  153
  154 U     153  155   86  154
  155 A     154  156   85  155
  156 G     155  157   84  156
  157 G     156  158   83  157
  158 G     157  159   82  158
  159 G     158  160   81  159
  160 U     159  161    0  160
  161 U     160  162    0  161
  162 G     161  163    0  162
  163 A     162  164    0  163
  164 A     163  165    0  164
  165 A     164  166    0  165
  166 C     165  167    0  166
  167 G     166  168   68  167
  168 C     167  169   67  168
  169 A     168  170   66  169
  170 G     169  171   65  170
  171 U     170  172    0  171
  172 C     171  173   61  172
  173 C     172  174   60  173
  174 U     173  175   30  174
  175 U     174  176   29  175
  176 C     175  177   28  176
  177 C     176  178   27  177
  178 C     177  179   26  178
  179 G     178  180   18  179
  180 G     179  181   17  180
  181 A     180  182   16  181
  182 G     181  183   15  182
  183 C     182  184    0  183
  184 A     183  185    0  184
  185 A     184  186    0  185
  186 G     185  187  226  186
  187 U     186  188  225  187
  188 A     187  189    0  188
  189 G     188  190  221  189
  190 G     189  191  220  190
  191 G     190  192  219  191
  192 G     191  193  218  192
  193 G     192  194  217  193
  194 G     193  195  216  194
  195 U     194  196  215  195
  196 C     195  197  214  196
  197 A     196  198    0  197
  198 A     197  199    0  198
  199 U     198  200  213  199
  200 G     199  201  212  200
  201 A     200  202  211  201
  202 G     201  203  210  202
  203 A     202  204    0  203
  204 A     203  205    0  204
  205 U     204  206   24  205
  206 G     205  207   23  206
  207 A     206  208   22  207
  208 U     207  209   21  208
  209 C     208  210   20  209
  210 U     209  211  202  210
  211 G     210  212  201  211
  212 A     211  213  200  212
  213 A     212  214  199  213
  214 G     213  215  196  214
  215 A     214  216  195  215
  216 C     215  217  194  216
  217 C     216  218  193  217
  218 U     217  219  192  218
  219 C     218  220  191  219
  220 C     219  221  190  220
  221 C     220  222  189  221
  222 U     221  223    0  222
  223 U     222  224    0  223
  224 G     223  225    0  224
  225 A     224  226  187  225
  226 C     225  227  186  226
  227 G     226  228    0  227
  228 C     227  229    0  228
  229 A     228  230    0  229
  230 U     229  231    0  230
  231 A     230  232    0  231
  232 G     231  233    0  232
  233 U     232  234    0  233
  234 C     233  235    0  234
  235 G     234  236    0  235
  236 A     235  237    0  236
  237 A     236  238    0  237
  238 U     237  239    0  238
  239 C     238  240    0  239
  240 C     239  241    0  240
  241 C     240  242    0  241
  242 C     241  243    0  242
  243 C     242  244    0  243
  244 A     243  245    0  244
  245 A     244  246    0  245
  246 A     245  247    0  246
  247 U     246  248    0  247
  248 a     247  249    0  248
  249 c     248  250    0  249
  250 a     249  251    0  250
  251 g     250  252    0  251
  252 a     251  253    0  252
  253 a     252  254    0  253
  254 g     253  255   14  254
  255 c     254  256   13  255
  256 g     255  257   12  256
  257 g     256  258   11  257
  258 g     257  259   10  258
  259 c     258  260    8  259
  260 u     259  261    7  260
  261 u     260  262    6  261
  262 a     261    0    0  262
"""
        self.not_enough_pseudoknot_tokens_path = TEST_DIR_PATH

        self.record_fd, self.record_path = tempfile.mkstemp()
        open(self.record_path, 'w').write(self.record_str)
        self.record = RNAStrandRecord(self.record_path)
        
        self.ref_defline = '> File ASE_00001.ct| RNA SSTRAND database, RNase \
P RNA, Acidianus ambivalens, strain DSM 3772| External source: RNase P Databas\
e, file name: A.ambivalens.ct2|'
        self.ref_sequence = 'GAGGAAAGUCCCGCCUCCAGAUCAAGGGAAGUCCCGCGAGG\
GACAAGGGUAGUACCCUUGGCAACUGCACAGAAAACUUACCCCUAAAUAUUCAAUGAGGAUUUGAUUCGACUCUUACCU\
UGGCGACAAGGUAAGAUAGAUGAAGAGAAUAUUUAGGGGUUGAAACGCAGUCCUUCCCGGAGCAAGUAGGGGGGUCAAU\
GAGAAUGAUCUGAAGACCUCCCUUGACGCAUAGUCGAAUCCCCCAAAUACAGAAGCGGGCUUA'
        self.ref_structure = ViennaStructure(\
                         '.....(((.(((((((((.[[[[[.((((((((((....)))))(((((..'+\
                         '...)))))((...((((............(((((((((((...........'+\
                         '.....((...((((((((((....))))))))))..))........)))))'+\
                         ')))))).......)))).)))))))))))...{{.{{{{{{{{..{{{{..'+\
                         ']]]]]}}}}}}}}}}}}...}}...........................))'+\
                         ')))))).')
        self.ref_vienna = self.ref_defline + '\n' + \
                          self.ref_sequence + \
                          '\n' + str(self.ref_structure) + '\n'
        
        self.ref_pairs = BasePairs([(6, 261), (7, 260), (8, 259), (10, 258),
            (11, 257), (12, 256), (13, 255), (14, 254), (15, 182), (16, 181),
            (17, 180), (18, 179), (20, 209), (21, 208), (22, 207), (23, 206),
            (24, 205), (26, 178), (27, 177), (28, 176), (29, 175), (30, 174),
            (31, 44), (32, 43), (33, 42), (34, 41), (35, 40), (45, 59),
            (46, 58), (47, 57), (48, 56), (49, 55), (60, 173), (61, 172),
            (65, 170), (66, 169), (67, 168), (68, 167), (81, 159), (82, 158),
            (83, 157), (84, 156), (85, 155), (86, 154), (87, 153), (88, 152),
            (89, 151), (90, 150), (91, 149), (108, 140), (109, 139),
            (113, 136), (114, 135), (115, 134), (116, 133), (117, 132),
            (118, 131), (119, 130), (120, 129), (121, 128), (122, 127),
            (186, 226), (187, 225), (189, 221), (190, 220), (191, 219),
            (192, 218), (193, 217), (194, 216), (195, 215), (196, 214),
            (199, 213), (200, 212), (201, 211), (202, 210)])

        # Scroll to the end of self.ref_not_en_tokens_pairs to see what's
        # happening to that list afterwards.
        self.ref_not_en_tokens_pairs = BasePairs([(2, 2918), (3, 2917),
            (4, 2916), (5, 2915), (11, 530), (12, 529), (13, 528), (14, 527),
            (15, 526), (16, 525), (17, 524), (18, 523), (19, 522), (20, 521),
            (26, 515), (27, 479), (28, 478), (29, 452), (31, 450), (32, 449),
            (33, 447), (34, 446), (35, 445), (36, 444), (37, 443), (38, 442),
            (39, 440), (40, 439), (41, 148), (42, 147), (44, 114), (46, 112),
            (47, 113), (48, 111), (49, 110), (50, 109), (51, 108), (52, 65),
            (53, 64), (54, 63), (56, 88), (57, 87), (59, 85), (60, 84),
            (61, 83), (62, 69), (71, 104), (72, 103), (73, 102), (74, 101),
            (75, 100), (76, 99), (77, 98), (80, 92), (81, 91), (82, 90),
            (115, 123), (116, 122), (117, 121), (130, 146), (131, 145),
            (132, 144), (133, 143), (134, 142), (135, 141), (151, 184),
            (152, 183), (153, 181), (154, 180), (155, 179), (156, 178),
            (157, 177), (163, 170), (164, 169), (167, 218), (191, 202),
            (193, 201), (194, 200), (203, 435), (204, 436), (205, 232),
            (206, 231), (207, 230), (208, 229), (209, 228), (215, 223),
            (216, 222), (235, 433), (236, 432), (237, 431), (238, 430),
            (243, 266), (244, 265), (245, 264), (248, 259), (249, 258),
            (250, 257), (251, 256), (272, 376), (273, 375), (274, 373),
            (275, 372), (276, 371), (277, 370), (278, 369), (279, 368),
            (286, 364), (287, 363), (288, 362), (289, 361), (296, 353),
            (297, 352), (298, 351), (299, 350), (300, 349), (301, 348),
            (302, 347), (306, 322), (307, 338), (308, 321), (309, 320),
            (310, 319), (311, 318), (323, 340), (324, 339), (325, 329),
            (331, 343), (332, 342), (333, 341), (381, 406), (382, 405),
            (383, 404), (384, 403), (385, 402), (386, 401), (387, 400),
            (388, 399), (389, 398), (390, 397), (393, 416), (411, 427),
            (412, 426), (413, 425), (414, 424), (415, 423), (417, 2448),
            (418, 2447), (419, 2446), (420, 2445), (421, 2444), (422, 2443),
            (462, 476), (464, 475), (465, 474), (466, 473), (467, 472),
            (480, 484), (489, 501), (490, 500), (491, 499), (492, 498),
            (534, 2081), (536, 2058), (538, 616), (539, 615), (540, 614),
            (541, 613), (542, 612), (543, 611), (544, 610), (545, 609), 
                (546, 608), (547, 607), (552, 605), (553, 604), (554, 601),
                (555, 600), (556, 599), (557, 598), (558, 597), (559, 596),
                (560, 595), (561, 594), (562, 593), (566, 590), (571, 584),
                (572, 583), (573, 582), (574, 581), (575, 580), (603, 606),
                (619, 634), (620, 633), (621, 632), (622, 631), (628, 2069),
                (635, 1364), (636, 1363), (637, 1362), (638, 1361), (639, 1360),
                (640, 1359), (642, 1354), (644, 760), (645, 758), (646, 757),
                (647, 756), (648, 755), (649, 754), (650, 753), (651, 752),
                (652, 751), (653, 750), (654, 749), (655, 748), (656, 747),
                (657, 746), (660, 684), (661, 683), (662, 682), (663, 681),
                (666, 678), (667, 677), (668, 676), (669, 675), (688, 695),
                (689, 694), (698, 726), (701, 725), (702, 724), (703, 723),
                (704, 722), (708, 718), (709, 717), (710, 716), (711, 715),
                (727, 742), (728, 741), (729, 740), (730, 739), (731, 738),
                (761, 901), (762, 900), (764, 898), (765, 895), (766, 894),
                (768, 891), (769, 890), (770, 889), (771, 888), (772, 887),
                (773, 886), (777, 867), (778, 866), (779, 865), (780, 864),
                (781, 863), (782, 862), (783, 861), (784, 860), (785, 859),
                (786, 858), (787, 857), (788, 855), (790, 822), (791, 821),
                (792, 820), (793, 819), (797, 814), (798, 813), (799, 812),
                (800, 811), (801, 810), (802, 809), (803, 808), (827, 852),
                (828, 851), (829, 850), (830, 849), (832, 848), (833, 847),
                (835, 846), (836, 845), (839, 2054), (869, 879), (870, 878),
                (871, 877), (904, 1299), (905, 1298), (906, 1297), (907, 1296),
                (908, 1295), (909, 1294), (910, 1293), (911, 1292), (912, 1070),
                (913, 1069), (914, 927), (915, 926), (916, 925), (917, 924),
                (918, 923), (929, 1040), (930, 1038), (931, 1037), (932, 1036),
                (933, 1035), (934, 1034), (935, 1033), (939, 1025), (940, 1024),
                (941, 1023), (943, 1022), (944, 1021), (945, 1020), (946, 1019),
                (947, 1018), (948, 1017), (949, 1016), (950, 1015), (951, 1014),
                (955, 1010), (956, 1008), (957, 1007), (963, 1003), (964, 1002),
                (965, 1001), (966, 1000), (967, 999), (968, 998), (969, 997),
                (970, 996), (971, 995), (972, 994), (974, 993), (975, 992),
                (976, 991), (977, 990), (978, 989), (979, 988), (980, 987),
                (981, 986), (1013, 2301), (1044, 1068), (1045, 1067),
                (1046, 1066), (1047, 1065), (1048, 1064), (1049, 1063),
                (1050, 1062), (1051, 1061), (1052, 1060), (1053, 1059),
                (1071, 1087), (1073, 1084), (1074, 1083), (1075, 1082),
                (1088, 1266), (1089, 1265), (1090, 1264), (1091, 1263),
                (1092, 1262), (1093, 1261), (1094, 1260), (1098, 1256),
                (1099, 1255), (1100, 1254), (1101, 1240), (1102, 1239),
                (1103, 1238), (1107, 1253), (1108, 1246), (1109, 1252),
                (1110, 1251), (1111, 1250), (1112, 1249), (1113, 1248),
                (1114, 1247), (1115, 1245), (1117, 1243), (1120, 1242),
                (1121, 1241), (1131, 1229), (1132, 1228), (1133, 1227),
                (1134, 1226), (1135, 1225), (1137, 1224), (1138, 1223),
                (1139, 1222), (1140, 1221), (1141, 1220), (1142, 1219),
                (1143, 1218), (1144, 1217), (1145, 1216), (1146, 1215),
                (1152, 2785), (1154, 1211), (1155, 1210), (1156, 1209),
                (1157, 1208), (1158, 1207), (1160, 1184), (1161, 1183),
                (1162, 1182), (1165, 1179), (1166, 1178), (1167, 1177),
                (1185, 1189), (1190, 1205), (1194, 1203), (1195, 1202),
                (1267, 1289), (1268, 1288), (1269, 1285), (1270, 1284),
                (1271, 1283), (1272, 1282), (1273, 1281), (1300, 1353),
                (1302, 1350), (1303, 1349), (1304, 1348), (1305, 1347),
                (1306, 1346), (1307, 1345), (1308, 1344), (1309, 1343),
                (1310, 1342), (1311, 1341), (1318, 1337), (1319, 1336),
                (1320, 1335), (1321, 1334), (1322, 1333), (1323, 1332),
                (1324, 1331), (1325, 1330), (1365, 2057), (1366, 2056),
                (1367, 2055), (1370, 2053), (1371, 2052), (1372, 2051),
                (1373, 2050), (1376, 1682), (1381, 1399), (1382, 1398),
                (1383, 1397), (1384, 1396), (1385, 1395), (1386, 1394),
                (1400, 1720), (1401, 1719), (1402, 1718), (1403, 1717),
                (1406, 1716), (1408, 1699), (1409, 1698), (1410, 1697),
                (1411, 1696), (1414, 1679), (1415, 1678), (1417, 1444),
                (1418, 1677), (1419, 1443), (1420, 1442), (1421, 1441),
                (1422, 1440), (1423, 1439), (1424, 1438), (1428, 1436),
                (1429, 1435), (1447, 1676), (1448, 1512), (1449, 1511),
                (1450, 1675), (1451, 1674), (1452, 1673), (1454, 1489),
                (1455, 1488), (1459, 1482), (1460, 1481), (1461, 1480),
                (1462, 1479), (1463, 1478), (1464, 1477), (1465, 1475),
                (1466, 1474), (1467, 1473), (1493, 1510), (1494, 1509),
                (1495, 1508), (1496, 1507), (1497, 1506), (1513, 1671),
                (1514, 1670), (1515, 1669), (1516, 1668), (1517, 1667),
                (1518, 1666), (1519, 1665), (1520, 1664), (1528, 1661),
                (1529, 1660), (1530, 1659), (1534, 1649), (1535, 1648),
                (1536, 1647), (1537, 1646), (1538, 1645), (1539, 1644),
                (1540, 1643), (1541, 1642), (1542, 1641), (1543, 1640),
                (1545, 1638), (1546, 1637), (1547, 1636), (1548, 1635),
                (1549, 1634), (1550, 1569), (1551, 1568), (1552, 1567),
                (1553, 1566), (1554, 1565), (1555, 1564), (1556, 1563),
                (1561, 2737), (1573, 1621), (1574, 1620), (1575, 1619),
                (1576, 1618), (1577, 1617), (1582, 1611), (1583, 1610),
                (1584, 1609), (1585, 1608), (1586, 1607), (1587, 1606),
                (1588, 1605), (1589, 1604), (1590, 1603), (1591, 1601),
                (1592, 1600), (1593, 1599), (1626, 1633), (1627, 1632),
                (1683, 1721), (1685, 1694), (1686, 1693), (1701, 1715),
                (1702, 1714), (1703, 1713), (1704, 1712), (1724, 2049),
                (1725, 2048), (1726, 2047), (1727, 2046), (1733, 2044),
                (1734, 2043), (1735, 2042), (1736, 2041), (1737, 2040),
                (1738, 2039), (1739, 2038), (1740, 2037), (1742, 2035),
                (1743, 2034), (1744, 2033), (1745, 2031), (1751, 2030),
                (1759, 1783), (1760, 1782), (1761, 1781), (1762, 1780),
                (1763, 1779), (1767, 1774), (1768, 1773), (1784, 1806),
                (1785, 1805), (1786, 1804), (1787, 1803), (1788, 1802),
                (1789, 1801), (1790, 1800), (1791, 1799), (1792, 1798),
                (1793, 1797), (1807, 1811), (1819, 2028), (1820, 2027),
                (1821, 2026), (1822, 2025), (1823, 2024), (1824, 2023),
                (1825, 2022), (1826, 2020), (1827, 2019), (1830, 1844),
                (1831, 1843), (1832, 1842), (1833, 1840), (1837, 2620),
                (1838, 2641), (1847, 1882), (1848, 1881), (1849, 1880),
                (1850, 1879), (1851, 1878), (1852, 1877), (1853, 1876),
                (1854, 1875), (1855, 1872), (1858, 1869), (1859, 1868),
                (1860, 1867), (1861, 1866), (1885, 2015), (1886, 2014),
                (1887, 2013), (1888, 2012), (1889, 2010), (1890, 1945),
                (1891, 1944), (1892, 1943), (1895, 1942), (1896, 1941),
                (1897, 1938), (1898, 1937), (1899, 1936), (1900, 1935),
                (1901, 1934), (1904, 1933), (1905, 1932), (1906, 1931),
                (1911, 1926), (1912, 1925), (1913, 1924), (1914, 1923),
                (1915, 1922), (1946, 1964), (1947, 1963), (1948, 1962),
                (1949, 1961), (1950, 1960), (1951, 1959), (1965, 1969),
                (1972, 2008), (1973, 2007), (1974, 2004), (1975, 2002),
                (1976, 1983), (1985, 2001), (1986, 2000), (1987, 1999),
                (1988, 1998), (1989, 1997), (1990, 1996), (1995, 2590),
                (2059, 2075), (2060, 2074), (2062, 2082), (2063, 2080),
                (2064, 2079), (2065, 2078), (2066, 2077), (2067, 2076),
                (2068, 2073), (2083, 2659), (2084, 2658), (2085, 2657),
                (2086, 2656), (2087, 2655), (2088, 2654), (2089, 2653),
                (2090, 2652), (2092, 2651), (2093, 2650), (2094, 2649),
                (2096, 2646), (2097, 2645), (2104, 2480), (2105, 2479),
                (2106, 2478), (2107, 2477), (2108, 2466), (2109, 2476),
                (2110, 2475), (2111, 2472), (2112, 2471), (2113, 2470),
                (2114, 2469), (2115, 2468), (2117, 2275), (2118, 2274),
                (2119, 2273), (2120, 2272), (2124, 2267), (2125, 2266),
                (2126, 2265), (2127, 2264), (2128, 2263), (2129, 2262),
                (2130, 2261), (2133, 2240), (2134, 2239), (2135, 2238),
                (2136, 2237), (2137, 2236), (2138, 2235), (2139, 2234),
                (2140, 2233), (2141, 2232), (2142, 2231), (2143, 2230),
                (2144, 2229), (2145, 2228), (2146, 2227), (2147, 2226),
                (2148, 2225), (2149, 2224), (2150, 2223), (2152, 2213),
                (2153, 2214), (2157, 2216), (2160, 2222), (2161, 2221),
                (2162, 2220), (2163, 2219), (2164, 2218), (2167, 2205),
                (2168, 2204), (2169, 2203), (2178, 2200), (2179, 2199),
                (2180, 2198), (2182, 2195), (2183, 2194), (2184, 2193),
                (2185, 2192), (2186, 2191), (2244, 2255), (2245, 2254),
                (2246, 2253), (2247, 2252), (2278, 2291), (2279, 2289),
                (2280, 2288), (2281, 2287), (2292, 2314), (2293, 2313),
                (2294, 2312), (2295, 2311), (2296, 2310), (2297, 2309),
                (2315, 2463), (2316, 2358), (2317, 2420), (2318, 2419),
                (2320, 2377), (2322, 2376), (2323, 2375), (2324, 2374),
                (2325, 2373), (2326, 2372), (2327, 2371), (2328, 2370),
                (2329, 2355), (2330, 2354), (2332, 2350), (2333, 2349),
                (2334, 2348), (2335, 2347), (2336, 2346), (2337, 2345),
                (2356, 2365), (2357, 2364), (2359, 2425), (2361, 2423),
                (2362, 2422), (2363, 2421), (2378, 2407), (2380, 2406),
                (2381, 2405), (2382, 2404), (2383, 2403), (2384, 2402),
                (2385, 2401), (2386, 2400), (2387, 2399), (2388, 2398),
                (2408, 2417), (2409, 2416), (2410, 2415), (2428, 2460),
                (2429, 2459), (2431, 2458), (2432, 2457), (2433, 2456),
                (2434, 2455), (2435, 2454), (2437, 2453), (2438, 2452),
                (2439, 2451), (2440, 2450), (2441, 2449), (2481, 2485),
                (2484, 2535), (2488, 2532), (2489, 2530), (2490, 2529),
                (2491, 2528), (2493, 2527), (2494, 2524), (2495, 2523),
                (2496, 2522), (2497, 2521), (2498, 2520), (2499, 2519),
                (2500, 2518), (2504, 2514), (2505, 2513), (2540, 2617),
                (2541, 2616), (2542, 2614), (2543, 2613), (2544, 2612),
                (2545, 2611), (2546, 2608), (2547, 2605), (2548, 2604),
                (2549, 2603), (2550, 2602), (2551, 2601), (2554, 2579),
                (2555, 2578), (2557, 2574), (2558, 2573), (2559, 2572),
                (2560, 2571), (2561, 2570), (2562, 2569), (2581, 2595),
                (2582, 2594), (2583, 2593), (2584, 2592), (2585, 2591),
                (2621, 2642), (2622, 2640), (2623, 2639), (2624, 2638),
                (2625, 2637), (2626, 2635), (2627, 2634), (2628, 2633),
                (2660, 2811), (2665, 2822), (2666, 2821), (2667, 2820),
                (2668, 2819), (2669, 2818), (2670, 2817), (2671, 2816),
                (2672, 2815), (2673, 2812), (2674, 2809), (2675, 2808),
                (2676, 2807), (2677, 2806), (2678, 2805), (2681, 2711),
                (2682, 2710), (2683, 2709), (2684, 2708), (2685, 2707),
                (2686, 2706), (2687, 2705), (2688, 2704), (2712, 2766),
                (2713, 2765), (2714, 2764), (2715, 2763), (2716, 2762),
                (2720, 2760), (2721, 2759), (2722, 2758), (2727, 2753),
                (2728, 2752), (2729, 2751), (2730, 2750), (2731, 2749),
                (2732, 2745), (2733, 2744), (2734, 2743), (2735, 2742),
                (2736, 2741), (2755, 2895), (2768, 2804), (2769, 2803),
                (2770, 2802), (2771, 2801), (2776, 2796), (2777, 2795),
                (2778, 2794), (2779, 2793), (2780, 2792), (2826, 2912),
                (2827, 2911), (2828, 2910), (2829, 2909), (2830, 2908),
                (2831, 2847), (2832, 2846), (2833, 2845), (2834, 2844),
                (2835, 2843), (2837, 2842), (2852, 2904), (2853, 2903),
                (2854, 2902), (2855, 2900), (2856, 2899), (2857, 2898),
                (2858, 2897), (2859, 2896), (2860, 2894), (2861, 2893),
                (2862, 2892), (2865, 2888), (2868, 2887), (2869, 2886),
                (2870, 2885), (2871, 2884), (2872, 2883), (2875, 2880)])
        
        self.ref_not_en_tokens_pairs = [(pair[0]+1, pair[1]+1) for pair in
            self.ref_not_en_tokens_pairs]

    def tearDown(self):
        os.close(self.record_fd)
        os.unlink(self.record_path)

class ComposeSsFromCmfinderMotivesTests(TestCase):
    """tests for compose_ss_from_cmfinder_motives
    """

    def setUp(self):
        motif_1 = 'test_data/output/cmfinder/2WDM_W/seq.fasta.1.motif.h1_1'
        motif_2 = 'test_data/output/cmfinder/2WDM_W/seq.fasta.1.motif.h1_1.h1_2'
        motif_3 = 'test_data/output/cmfinder/2WDM_W/seq.fasta.1.motif.h1_2'
        motif_4 = 'test_data/output/cmfinder/2WDM_W/seq.fasta.1.motif.h1_3'
        motif_5 = 'test_data/output/cmfinder/2WDM_W/seq.fasta.1.motif.h2_1'
        motif_6 = 'test_data/output/cmfinder/2WDM_W/seq.fasta.1.motif.h2_1.h2_2'
        motif_7 = 'test_data/output/cmfinder/2WDM_W/seq.fasta.1.motif.h2_2'
        
        self.motives = [open(x) for x in (motif_1, motif_2, motif_3, motif_4,
                                          motif_5, motif_6, motif_7)]

    def test_2wdm_w(self):
        """compose_ss_from_cmfinder_motives should work correctly on 2WDM_W
        """
        input_seq = 'GCCCGGAUAGCUCAGUCGGUAGAGCAGGGGAUUGAAAAUCCCCGUGUC'+\
                    'CUUGGUUCGAUUCCGAGUCCGGGCACCA'
        ss = compose_ss_from_cmfinder_motives(self.motives, '2WDM_W',
                                              input_seq)
        
        self.assertEqual(len(ss), len(input_seq))
        self.assertEqual(ss,
                         '.((((((...(((........)))(((((((((..'+\
                         '.))))))).))..(((((.......))))))))))).....')
        
    def test_3FIH_Y(self):
        """compose_ss_from_cmfinder_motives should work on 3FIH_Y
        """
        input_seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUG'+\
                    'UGUUCGAUCCACAGAAUUCGCACCA'
        ss = compose_ss_from_cmfinder_motives(self.motives, '3FIH_Y',
                                              input_seq)
        
        self.assertEqual(len(ss), len(input_seq))
        self.assertEqual(ss,
                         '.((((((...(((........)))(.(((((.(...).)))))..).'+\
                         '.(((((.......))))))))))).....')
        
    def test_corynebcaterium(self):
        """compose_ss_from_cmfinder_motives should work on Corynebacterium seqs
        """
        motif_1 = \
            'test_data/output/cmfinder/corynebacterium/seq.fasta.1.motif.h1_1'
        motif_2 = \
            'test_data/output/cmfinder/corynebacterium/seq.fasta.1.motif.h2_1'

        motives = [open(x) for x in (motif_1, motif_2)]
        
        input_seq = 'CGAGUUGCAACGGGAGAGCUUGCUUUCCGCGCCGGGCGCAUUGUGCUUCCGGCCA'+\
                    'GGAAAGUAGCGCCGUAGGAGCAAACCCUCCCCGGGAAUCUCUCAGGCCCAGUACC'+\
                    'GCUGCAA'
        ss = compose_ss_from_cmfinder_motives(motives, 'Corynebacterium_krop.1',
                                              input_seq)
        self.assertEqual(len(ss), len(input_seq))
        self.assertEqual(ss,
                         '...............................................'+\
                         '..................(((...((((....(((.....)))...)'+\
                         ')))..)))...............')

    def test_glmS(self):
        """compose_ss_from_cmfinder_motives should work on glmS seqs
        """
        motif_1 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h1_1'
        motif_2 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h1_1.h2_1'
        motif_3 = \
            'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h1_1.h2_1.h2_3'
        motif_4 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h1_2'
        motif_5 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h1_2.h2_2'
        motif_6 = \
            'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h1_2.h2_2.h1_3'
        motif_7 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h1_3'
        motif_8 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h2_1'
        motif_9 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h2_2'
        motif_10 = 'test_data/output/cmfinder/glmS/seq.fasta.1.motif.h2_3'
        
        motives = [open(x) for x in (motif_1, motif_2, motif_3, motif_4,
                                     motif_5, motif_6, motif_7, motif_8,
                                     motif_9, motif_10)]
        
        input_seq = 'CGTTGTATGTTGCTTGCGCATATTTTGTAGTAGTTGTTTGAAAATGTTGA'+\
                    'TCAGTCTGAATGAAACAAAAAAAGTTTGGAGATTTGTTGGGAATGATTGA'+\
                    'CGTGCAGACAACAATTAGGTATGATTAGCCTGTCGAAGATAGAAAGGAGG'+\
                    'ATCGCAACAGAGGGAGACCTCTATGACAGTTGAAGCGCCAGGACTATTCC'+\
                    'ACGGACGGTAACGGAATAGTTGACGAGGAGAGGGTTTATCGAAGGTTCGG'+\
                    'CGGATGCCCTCCGGTTGCACATGACAGCCGCAAGCTTTTGTAAAAAACAG'+\
                    'AGGGGCGACCTTCTGGACAAAGGCAAAAGCGTTCATGCTAACAAAAAAGA'+\
                    'GATAAAGAGTGGGGCAAG'
        
        ss = compose_ss_from_cmfinder_motives(motives,
                                              'AP001507.1/288594-288961',
                                              input_seq)
        self.assertEqual(len(ss), len(input_seq))
        self.assertEqual(ss,
                         '...................................................'+\
                         '...................................................'+\
                         '...................................................'+\
                         '..................................................('+\
                         '((.((......)).)))........((((((...(((......))).....'+\
                         '))))))((((((.......))))))....(((((((......(((((((..'+\
                         '..))))))).......)))))))............................'+\
                         '...........')

if __name__ == "__main__":
    main()
