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

import os
ABSPATH = os.path.split(os.path.abspath(__file__))[0]

from cogent import LoadSeqs
from cogent.util.unit_test import TestCase, main
from cogent.struct.pairs_util import correlation_coefficient, extract_seqs
from cogent.core.sequence import RnaSequence, ModelSequence
from cogent.core.moltype import RNA

from .secstruc import BasePairs
from .utils import mcc, selectivity, sensitivity, get_counts, \
make_sequence_alignment, extract_query_seq_pred, ExtractQuerySeqPredError,\
get_all_pairs, get_bps_for_aligned_seq
from .utils import CannotMakeSequenceAlignmentError as CMSAE

class GardnerMetricsTest(TestCase):
    """Tests for the metrics from Gardner & Giegerich 2004"""
    
    def setUp(self):
        """setUp: setup method for all tests"""
        self.true = BasePairs([(0,40),(1,39),(2,38),(3,37),(10,20),\
            (11,19),(12,18),(13,17),(26,33),(27,32)])
        self.predicted = BasePairs([(0,40),(1,39),(2,38),(3,37),(4,36),\
            (5,35),(10,22),(11,20),(14,29),(15,28)])
        self.seq = ['>seq1\n','agguugaaggggauccgauccacuccccggcuggucaaccu']

    def test_conflicts(self):
        """all metrics should raise error when conflicts in one of the structs
        """
        ref = BasePairs([(1,6),(2,5),(3,10),(7,None),(None,None),(5,2),(1,12)])
        pred = BasePairs([(6,1),(10,11),(3,12)])
        
        self.assertRaises(ValueError, sensitivity, ref, pred)
        self.assertRaises(ValueError, sensitivity, pred, ref)
        self.assertRaises(ValueError, selectivity, ref, pred)
        self.assertRaises(ValueError, selectivity, pred, ref)
        self.assertRaises(ValueError, mcc, ref, pred, self.seq)
        self.assertRaises(ValueError, mcc, pred, ref, self.seq)

    def test_get_all_pairs(self):
        """get_all_pairs: should return the number of possible pairs"""
        seq = RnaSequence('UCAG-NACGU')
        seq2 = RnaSequence('UAAG-CACGC')
        self.assertEqual(get_all_pairs([seq], min_dist=4), 6)
        self.assertEqual(get_all_pairs([seq2], min_dist=4), 4)
        # when given multiple sequences, should average over all of them
        self.assertEqual(get_all_pairs([seq,seq2], min_dist=4), 5)
        # different min distance
        self.assertEqual(get_all_pairs([seq], min_dist=2),10)
        # error on invalid minimum distance
        self.assertRaises(ValueError, get_all_pairs, [seq], min_dist=-2)

    def test_get_counts(self):
        """get_counts: should work with all parameters"""
        seq = RnaSequence('UCAG-NAUGU')
        p = BasePairs([(1,8),(2,7)])
        p2 = BasePairs([(1,8),(2,6),(3,6),(4,9),])
        exp = {'TP':1,'TN':0, 'FN':1,'FP':3,\
            'FP_INCONS':0, 'FP_CONTRA':0, 'FP_COMP':0}
        self.assertEqual(get_counts(p, p2, False), exp)
        exp = {'TP':1,'TN':0, 'FN':1,'FP':3,\
            'FP_INCONS':1, 'FP_CONTRA':1, 'FP_COMP':1}
        self.assertEqual(get_counts(p, p2, split_fp=True), exp)
        seq = RnaSequence('UCAG-NACGU')
        exp = {'TP':1,'TN':7, 'FN':1,'FP':3,\
            'FP_INCONS':1, 'FP_CONTRA':1, 'FP_COMP':1}
        self.assertEqual(get_counts(p, p2, split_fp=True,\
            sequences=[seq], min_dist=2), exp)
        # check against compare_ct.pm
        exp = {'TP':4,'TN':266, 'FN':6,'FP':6,\
            'FP_INCONS':2, 'FP_CONTRA':2, 'FP_COMP':2}
        seq = 'agguugaaggggauccgauccacuccccggcuggucaaccu'.upper()
        self.assertEqual(get_counts(self.true, self.predicted, split_fp=True,\
            sequences=[seq], min_dist=4), exp)

    def test_extract_seqs(self):
        """extract_seqs: should handle different input formats"""
        s1 = ">seq1\nACGUAGC\n>seq2\nGGUAGCG"
        s2 = [">seq1","ACGUAGC",">seq2","GGUAGCG"]
        s3 = ['ACGUAGC','GGUAGCG']
        s4 = [RnaSequence('ACGUAGC'), RnaSequence('GGUAGCG')]
        m1 = ModelSequence('ACGUAGC', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        m2 = ModelSequence('GGUAGCG', Name='rna2',\
            Alphabet=RNA.Alphabets.DegenGapped)
        s5 = [m1, m2]
        f = extract_seqs
        self.assertEqual(f(s1), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s2), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s3), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s4), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s5), ['ACGUAGC', 'GGUAGCG'])

    def test_sensitivity(self):
        """sensitivity: check against compare_ct.pm"""
        sen = sensitivity(self.true,self.predicted)
        self.assertEqual(sen, 0.4)

    def test_sensitivity_general(self):
        """sensitivity: should work in general"""
        ref = BasePairs([(1,6),(2,5),(3,10)])
        pred = BasePairs([(6,1),(10,11),(3,12)])
        # one good prediction
        self.assertFloatEqual(sensitivity(ref, pred), 1/3)
        # over-prediction not penalized
        pred = BasePairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        self.assertFloatEqual(sensitivity(ref, pred), 1/3)

    def test_sensitivity_dupl(self):
        """sensitivity: should handle duplicates, pseudo, None"""
        ref = BasePairs([(1,6),(2,5),(3,10),(7,None),(None,None),(5,2),(4,9)])
        pred = BasePairs([(6,1),(10,11),(3,12)])
        self.assertFloatEqual(sensitivity(ref, pred), 0.25)

        pred = BasePairs([(6,1),(10,11),(3,12),(20,None),(None,None),(1,6)])
        self.assertFloatEqual(sensitivity(ref, pred), 0.25)

    def test_sensitivity_empty(self):
        """sensitivity: should work on emtpy BasePairs"""
        # both empty
        self.assertFloatEqual(sensitivity(BasePairs([]), BasePairs([])), 1)
        pred = BasePairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        # prediction emtpy
        self.assertFloatEqual(sensitivity(BasePairs([]), pred), 0)
        # reference empty
        self.assertFloatEqual(sensitivity(pred, BasePairs([])), 0)

    def test_selectivity(self):
        """selectivity: check against compare_ct.pm"""
        sel = selectivity(self.true,self.predicted)
        self.assertEqual(sel, 0.5)

    def test_selectivity_general(self):
        """selectivity: should work in general"""
        ref = BasePairs([(1,6),(2,5),(10,13)])
        pred = BasePairs([(6,1),(3,4),(10,12)])
        # one good prediction
        self.assertFloatEqual(selectivity(ref, pred), 0.5)
        # over-prediction not penalized
        pred = BasePairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        self.assertFloatEqual(selectivity(ref, pred), 0.25)

    def test_selectivity_dupl(self):
        """selectivity: duplicates and Nones shouldn't influence the calc.
        """
        ref = BasePairs([(1,6),(2,5),(10,13),(6,1),(7,None),(None,None)])
        pred = BasePairs([(6,1),(3,4),(10,12)])
        self.assertFloatEqual(selectivity(ref, pred), 0.5)

    def test_selectivity_empty(self):
        """selectivity: should handle empty reference/predicted structure"""
        # both empty
        self.assertFloatEqual(selectivity(BasePairs([]), BasePairs([])), 1)
        pred = BasePairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        # prediction emtpy
        self.assertFloatEqual(selectivity(BasePairs([]), pred), 0)
        # reference empty
        self.assertFloatEqual(selectivity(pred, BasePairs([])), 0)

    def test_correlation_coefficient(self):
        """correlation_coefficient: check against compare_ct.pm"""
        self.assertFloatEqual(correlation_coefficient(self.true,\
            self.predicted, seqs=self.seq, min_dist=4), 0.42906394)

    def test_cc_bad_pred(self):
        """correlation_coefficient: should give 0 when TP=0"""
        ref = BasePairs([(1,7),(2,5)])
        pred = BasePairs([(0,8)])
        seqs = ['CAUCGAUUG']
        self.assertEqual(correlation_coefficient(ref, pred, seqs=seqs), 0.0)

    def test_mcc(self):
        """mcc: check against compare_ct.pm"""
        res = mcc(self.true,self.predicted,self.seq, min_dist=4)
        self.assertFloatEqual(res, 0.42906394)

    def test_get_counts_pseudo(self):
        """get_counts: should work when pseudo in ref -> classification off"""
        # pairs that would normally be compatible, are now contradicting
        ref = BasePairs([(0,8),(1,7),(4,10)])
        pred = BasePairs([(0,8),(3,6),(4,10)])
        seq = 'GACUGUGUCAU'
        exp = {'TP':2,'TN':13-2-1, 'FN':1,'FP':1,\
            'FP_INCONS':0, 'FP_CONTRA':1, 'FP_COMP':0}
        self.assertEqual(get_counts(ref, pred, split_fp=True,\
            sequences=[seq], min_dist=4), exp)

    def test_get_all_pairs_2(self):
        """get_all_pairs: should return the number of possible pairs"""
        seq = RnaSequence('UCAG-NACGU')
        seq2 = RnaSequence('UAAG-CACGC')
        self.assertEqual(get_all_pairs([seq], min_dist=4), 6)
        self.assertEqual(get_all_pairs([seq2], min_dist=4), 4)
        # when given multiple sequences, should average over all of them
        self.assertEqual(get_all_pairs([seq,seq2], min_dist=4), 5)
        # different min distance
        self.assertEqual(get_all_pairs([seq], min_dist=2),10)
        # error on invalid minimum distance
        self.assertRaises(ValueError, get_all_pairs, [seq], min_dist=-2)

class MakeSequenceAlignmentTests(TestCase):
    """tests for make_sequence_alignment function from output_unified.py
    """
    
    def setUp(self):
        """Sets up environment for tests
        """
        self.random_seq = LoadSeqs(data=\
        '>seq0\nACUGCGCGGAUCGAUCGAUCGAUCGAUGCAUUUUACGAUCGCCA\n', aligned=False)
        
        self.rrna = LoadSeqs(data=RRNA, aligned=False)
        self.rrna_aln = LoadSeqs(data=REF_ALN)
        self.seq_db_path = os.path.join(ABSPATH, 'test_data', 'Rfam10_part.fasta')
        
    def test_weired_input(self):
        """make_sequence_alignment raises error when input is incorrect"""
        self.assertRaises(CMSAE, make_sequence_alignment, '>seq\nACUG', \
                          self.seq_db_path)
        self.assertRaises(CMSAE, make_sequence_alignment, [], self.seq_db_path)

    def test_random_seq(self):
        """make_sequence_alignment raises error when Paralign fails"""
        self.assertRaises(CMSAE, make_sequence_alignment, self.random_seq, \
                          self.seq_db_path)
        
    def test_rrna(self):
        """make_sequence_alignment works correctly on correct input"""
        alignment = make_sequence_alignment(self.rrna, self.seq_db_path,
                                            max_nr_seqs=20)
        self.assertEqual(alignment, self.rrna_aln)
        
    def test_seq_coll(self):
        """make_sequence_alignment works correctly for input sequence collection
        """
        fasta_coll = ">seq0\nUACGGCGGUCAUAGCGGGCGUUCCGAACCCGGU\
CGUUAAGCCGCCCAGCGCCGAUGGUACUGCCCUGUGGUGGGGUGGGAGAGUAGGACACCGCCGUAC\n\
>seq1\nUACGGCGGUCAUAGCGGGCCCACCCAACCCGGU\
CGUUAAGCCGCCCAGCGCCGAUGGUAUGCCCUGUGGUGGGGUGGGAGAGUAGGACACCGCCGUAC"
        coll = LoadSeqs(data=fasta_coll, aligned=False)
        aln = make_sequence_alignment(coll, self.seq_db_path,
                                      max_nr_seqs=20)
        
        ref_aln = """>seq0  _F_ seq_.rfold
UACGGCGGUCAUAGCGGGCGUUCCGAACCCGGUCGUUAAGCCGCCCAGCGCCGAUGGUACUGCCCUGUGGUGGGGUGGGAGAGUAGGACACCGCCGUAC
>seq1  _F_ seq1.rfold
UACGGCGGUCAUAGCGGGCCCACCCAACCCGGUCGUUAAGCCGCCCAGCGCCGAUGGUA-UGCCCUGUGGUGGGGUGGGAGAGUAGGACACCGCCGUAC
"""

        self.assertEqual(aln, ref_aln)


RRNA = '>seq0\nUACGGCGGUCAUAGCGGGCGUUCCGAACCCGGU\n\
CGUUAAGCCGCCCAGCGCCGAUGGUACUGCCCUGUGGUGGGGUGGGAGAG\n\
UAGGACACCGCCGUAC\n'
        
REF_ALN = """>seq0  _F_ seq_.rfold
UACGGCGGUCAUAGC-GG-----------GC------GUUCCGAACCCGG
UCGUUAAGCCGCCCAGCGCCGAUGGUACUGCCCU-GUGGUGGGGUGGGAG
AGUAGGACAC-CG-CCGUAC
>RF00001_5S_rRNA_CP000088.1/2386815-2386930  _F_ RF____1_5S_rRNA_CP____88.1_2386815-238693_.rfold
UACGGCGGUCAUAGC-GGGCGGGGUCCACCCGGUCCCAUUCCGAACCCGG
UCGUUAAGCCGCCCAGCGCCGAUGGUACUGCCCU-GUGGUGGGGUGGGAG
AGUAGGACAC-CG-CCGUAC
>RF00001_5S_rRNA_D13617.1/2-117  _F_ RF____1_5S_rRNA_D13617.1_2-117.rfold
UGUGGUGGUUAUUGC-UGGAGGGUCACACCCGGUCUCUUUCCGAACCCGG
UCGUUAAGUCUCUUUGCGCUGAUGGUACUGCCCU-GUGGUGGGGUGGGAG
AGUAGGUCGC-UG-CCACGG
>RF00001_5S_rRNA_ABUA01000019.1/261-377  _F_ RF____1_5S_rRNA_ABUA_1____19.1_261-377.rfold
UUCGGUGGUCAUUGC-GGUUGGGGAACACCCGGUCCCAUUCCGAACCCGG
UAGUUAAGCCUUCCAGCGCCGAUGGUACUGCACUCGUGAGGGUGUGGGAG
AGUAGGACGC-CG-CCGAAC
>RF00001_5S_rRNA_ABUU01000012.1/99629-99745  _F_ RF____1_5S_rRNA_ABUU_1____12.1_99629-99745.rfold
UACGGCGGUCAUGGC-GAGAGGGAAACGCCCGGUCCCAUUCCGAACCCGG
AAGCUAAGCCUCUCAGCGCCGAUGGUACUGCCCUGGAGACGGGGUGGGAG
AGUAGGACAC-CG-CCGGAC
>RF00001_5S_rRNA_D13616.1/3-119  _F_ RF____1_5S_rRNA_D13616.1_3-119.rfold
UACGGCGGUUAUAGC-GGUAGGGAAACACCCGGUUACAUUCCGAACCCGG
AAGUUAAGCCUUCCAGCGCCGAUGGUACUGCACUGGGGACGGUGUGGGAG
AGUAGGACGC-CG-CCGGAA
>RF00001_5S_rRNA_ABUZ01000030.1/33258-33374  _F_ RF____1_5S_rRNA_ABUZ_1____3_.1_33258-33374.rfold
UCCGGUGGUUAUGG-CGGAGGGGAAACACCCGGUCCCAUUCCGAACCCGG
UAGUUAAGCCCUCCAGCGCCGAUGGUACUGCAUGGGAGACUGUGUGGGAG
AGUAGGUCAC-CG-CCGGAA
>RF00001_5S_rRNA_Z50076.1/4-120  _F_ RF____1_5S_rRNA_Z5__76.1_4-12_.rfold
GUCGGUGGCGACAGC-GGGGAGGGUCCGCCCGGUCCCAUACCGAACCCGG
AAGCUAAGCUCCCCAGCGCCGAUGGUACUGCACCCGGUGGGGUGUGGGAG
AGUAGGACAC-CG-CCGACA
>RF00001_5S_rRNA_ACGK01000003.1/5074-5189  _F_ RF____1_5S_rRNA_ACGK_1_____3.1_5_74-5189.rfold
GUCAGCGAUUUUGGC-AGGGGAGGUACACCCGGUCCCAUUCCGAACCCGG
AAGUUAAGCCCCCUAGCGCCGAUGGUACUGCGGA-GUUAGUCCGUGGGAG
AGCAGGACGU-CG-CUGACC
>RF00001_5S_rRNA_ABXR01000030.1/47727-47843  _F_ RF____1_5S_rRNA_ABXR_1____3_.1_47727-47843.rfold
CCGGGUGCCGAUACUGGGGCGGGAAACACCCGGUUCCAUUCCGAACCCGG
CCGUUAAGCCGCCCUGGGCCGAUGGUAGUAUGGGGGUAGCCCCAUGCGAG
AGUAGGUAG--CGCCC-GGG
>RF00001_5S_rRNA_AARF01000019.1/827-943  _F_ RF____1_5S_rRNA_AARF_1____19.1_827-943.rfold
UUUGGUGACGAUAGC-GGAAGGGAACCACGCGUACCCAUCCCGAACACGG
UCGUUAAGCCUUCCAGCGCCGAUGGUACUUGGACCCUAGGGUCCUGGGAG
ACUAGGACGU-UG-CCAAGC
>RF00001_5S_rRNA_ABZW01000010.1/314-430  _F_ RF____1_5S_rRNA_ABZW_1____1_.1_314-43_.rfold
UCCGGUGGCCAUAGU-GGAAGGGAAACACCCGGUCCCAUUCCGAACCCGG
UCGUUAAGCCUUCCAACGCUGAUGGUACUGCACGAGGGAUCGUGUGGGAG
AGUAAGACGC-UG-CCGGAA
>RF00001_5S_rRNA_U83911.1/5199-5316  _F_ RF____1_5S_rRNA_U83911.1_5199-5316.rfold
CACGGUGGCCAUAGC-GGAAGGGAAACACCCGGUCACAUUCCGAACCCGG
AAGUUAAGCCUUCCAGCGCCGAUGGUACUGCACCCGGGAGGGUGUGGGAG
AGUAGGUCGCUUG-CCGGGC
>RF00001_5S_rRNA_AE000657.1/1196967-1197083  _F_ RF____1_5S_rRNA_AE___657.1_1196967-1197_83.rfold
CCUGGUGGCCAUAGCGGGGG-GGAAACACCCGGUCCCAUUCCGAACCCGG
CAGUUAAGCCCCCCAGCGCCGAUGAUACUGUGCCGGCAGCGGCACGGGAA
AGUAGGUCGC-UG-CCAGGG
>RF00001_5S_rRNA_ABUH01000013.1/102410-102526  _F_ RF____1_5S_rRNA_ABUH_1____13.1_1_241_-1_2526.rfold
UUCGGCGGUUAUAGC-GUGAGGGAAACGCCCGGUCCCAUUCCGAACCCGG
AAGCUAAGCUUCACUGCGCCGAUGGUACUGCACCUGCGAGGGUGUGGGAG
AGUAGGACAC-CG-CCGGAC
>RF00001_5S_rRNA_ACJX01000067.1/5058-5173  _F_ RF____1_5S_rRNA_ACJX_1____67.1_5_58-5173.rfold
UAUGGUGGCCAUAGC-GGAGGGGAAACACCCGGAUCCAUGCCGAACCCGG
CAGUUAAGCCCUCCAGCGCCGAUGGUACUGCGGGGGGUA-CCCGUGGGAG
AGUAGGUCGC-UG-CCAUGU
>RF00001_5S_rRNA_X52302.1/2-118  _F_ RF____1_5S_rRNA_X523_2.1_2-118.rfold
GCUGGUGGCUAUGGC-GGAAGGGCCACACCCGAUCCCAUCCCGAACUCGG
UCGUUAAGACUUCCAGCGCCGAUGGUACUGCAGGGGCGACCCUGUGGGAG
AGUAGGACGC-CG-CCAGCG
>RF00001_5S_rRNA_Z50060.1/3-119  _F_ RF____1_5S_rRNA_Z5__6_.1_3-119.rfold
UUCGGUGGUUAUGGC-GGUAGGGAAACGCCCGGUCCCAUUCCGAACCCGG
AAGCUAAGCCUGCCAGCGCCGAUGGUACUGCACCCCCGCGGGUGUGGGAG
AGUAGGACAC-CG-CCGAAC
>RF00001_5S_rRNA_ABUD01000038.1/14311-14427  _F_ RF____1_5S_rRNA_ABUD_1____38.1_14311-14427.rfold
UUCGGUGGUCAUAGC-GAAGGGGAAACGCCCGGUUACAUUCCGAACCCGG
AAGCUAAGCCCUUCAGCGCCGAUGGUACUGCACUGGGAACGGUGUGGGAG
AGUAGGACAC-CG-CCGAAC
>RF00001_5S_rRNA_AP006618.1/2168794-2168910  _F_ RF____1_5S_rRNA_AP__6618.1_2168794-216891_.rfold
UACGGCGGUCAUAGC-GGUGGGGAAACGCCCGGUCCCAUUCCGAACCCGG
AAGCUAAGCCUGCCAGCGCCGAUGGUACUGCACUCGACAGGGUGUGGGAG
AGUAGGACAC-CG-CCGGAA
>RF00001_5S_rRNA_ACFH01000038.1/5466-5582  _F_ RF____1_5S_rRNA_ACFH_1____38.1_5466-5582.rfold
CUCGGUGGUCAUAGC-GGGGGAGCCACGCCCGGCCCCAUUCCGAACCCGG
AAGCUAAGACCCCCAGCGCCGAUGGUACUGCACCCGCCAGGGUGUGGGAG
AGUAGGACAC-CG-CCGAGC"""

class ExtractQuerySeqPredTests(TestCase):
    """tests for extract_query_seq_pred function from output_unified.py
    """
    
    def setUp(self):
        """setting up environment for tests
        """
        self.results = [
            {'vienna_sec_structs' : '>seq 0\nACUG\n....\n>seq 1\nACCC\n(..)\n',
             'program_name'       : 'method1', 'blablab' : 1,
             'sequence' : 'not important during this test',
             'all_energies' : []},
            {'vienna_sec_structs' : '>seq 0\nACUG\n(())\n>seq 1\nACCC\n.().\n',
             'program_name'       : 'stefan', 'blablab' : 2,
             'sequence' : 'not important during this test',
             'all_energies' : []},
            {'best_predicted_ss' : '....(((...)))....',
             'program_name' : 'test', 'blablabla' : 3,
             'sequence' : 'not important at all during this test',
             'all_energies' : [100.0], 'best_predicted_energy' : 100.0}
        ]
        
    def test_extract_query_seq_pred_ok(self):
        """extract_query_seq_pred should work correctly if correct args
        """
        first_result = extract_query_seq_pred(self.results[0], 'ACUG', 'seq 0')
        self.assertTrue(type(first_result) is dict)
        self.assertTrue('sequence' in first_result)
        self.assertTrue('program_name' in first_result)
        self.assertEqual(first_result.get('best_predicted_ss'), '....')
        
        second_result = extract_query_seq_pred(self.results[0], 'ACCC', 'seq 1')
        self.assertTrue(type(second_result) is dict)
        self.assertTrue('sequence' in second_result)
        self.assertTrue('program_name' in second_result)
        self.assertEqual(second_result.get('best_predicted_ss'), '(..)')
        
        third_result = extract_query_seq_pred(self.results[1], 'ACCC', 'seq 1')
        self.assertTrue(type(third_result) is dict)
        self.assertTrue('sequence' in third_result)
        self.assertTrue('program_name' in third_result)
        self.assertEqual(third_result.get('best_predicted_ss'), '.().')
        
    def test_extract_query_seq_pred_not_ok(self):
        """extract_query_seq_pred should work correctly if wrong args
        """
        self.assertRaises(ExtractQuerySeqPredError, extract_query_seq_pred,
                          self.results[0], 'ACCG', 'seq 0') # zła sekwencja
        
        self.assertRaises(ExtractQuerySeqPredError, extract_query_seq_pred,
                          self.results[0], 'ACUG', 'seq 1') # zła nazwa

    def test_extract_query_seq_pred_best_predicted_ss(self):
        """extract_query_seq_pred should work correctly for best_predicted_ss
        """
        first_result = extract_query_seq_pred(self.results[2], 'ACUG', 'seq 0')
        self.assertTrue(type(first_result) is dict)
        self.assertTrue('sequence' in first_result)
        self.assertTrue('program_name' in first_result)
        self.assertEqual(first_result.get('best_predicted_ss'),
                         '....(((...)))....')
        self.assertEqual(first_result.get('best_predicted_energy'),
                         100.0)

class GetBpsForAlignedSeqTests(TestCase):
    """tests for get_bps_for_aligned_seq function
    """
    
    def test_seq_simple(self):
        """get_bps_for_aligned_seq should work for simple case
        """
        aln_seq = "--U--------------A"
        pred = BasePairs([(2, 17)])
        
        result = get_bps_for_aligned_seq(aln_seq, pred)
        self.assertEqual(result, BasePairs([(0, 1)]))

    def test_seq_simple_offset(self):
        """get_bps_for_aligned_seq should work when first_index=1, part 1
        """
        aln_seq = "--U--------------A"
        pred = BasePairs([(2, 17)])
        
        result = get_bps_for_aligned_seq(aln_seq, pred, 1)
        self.assertEqual(result, BasePairs([]))

    def test_simple_offset_ok(self):
        """get_bps_for_aligned_seq should work when first_index=1, part 2
        """
        aln_seq = "-U--------------A-"
        pred = BasePairs([(2, 17)])
        
        result = get_bps_for_aligned_seq(aln_seq, pred, 1)
        self.assertEqual(result, BasePairs([(1, 2)]))
        
    def test_seq_conflict(self):
        """get_bps_for_aligned_seq should work for conflicted case
        """
        aln_seq = "--U--A-----------A"
        pred = BasePairs([(2, 5), (2, 17)])
        
        result = get_bps_for_aligned_seq(aln_seq, pred)
        self.assertEqual(result, BasePairs([(0, 1), (0, 2)]))
        
    def test_seq_conflict_offset(self):
        """get_bps_for_aligned_seq should conflict first_index=1, part 1
        """
        aln_seq = "--U--A-----------A"
        pred = BasePairs([(2, 5), (2, 17)])
        
        result = get_bps_for_aligned_seq(aln_seq, pred, 1)
        self.assertEqual(result, BasePairs([]))

    def test_seq_conflict_offset_ok(self):
        """get_bps_for_aligned_seq should conflict first_index=1, part 2
        """
        aln_seq = "-U--A-----------A-"
        pred = BasePairs([(2, 5), (2, 17)])
        
        result = get_bps_for_aligned_seq(aln_seq, pred, 1)
        self.assertEqual(result, BasePairs([(1,2), (1,3)]))
        
        result = get_bps_for_aligned_seq(aln_seq, pred, 2)
        self.assertEqual(result, BasePairs([]))
        
    def test_aligned_seq_skip(self):
        """get_bps_for_aligned_seq should work when some base pairs are skipped
        """
        aln_seq = "ACUAGCUG-----ACUGA"
        pred = BasePairs([(2, 10), (2, 17)])
        
        result = get_bps_for_aligned_seq(aln_seq, pred)
        self.assertEqual(result, BasePairs([(2,12)]))
        
    def test_many_gaps_seq(self):
        """get_bps_for_aligned_seq test on seq with many gaps
        """
        aln_seq = "ACUAGCUG-----ACU-A---------CGCGC---A"
        pred = BasePairs([(2, 10), (2, 17), (4, 36)])
        
        result_offset = get_bps_for_aligned_seq(aln_seq, pred, 1)
        self.assertEqual(result_offset, BasePairs([(4, 18)]))
        
        result_offset = get_bps_for_aligned_seq(aln_seq, pred, 2)
        self.assertEqual(result_offset, BasePairs([(2, 12)]))

if __name__ == "__main__":
    main()
    
