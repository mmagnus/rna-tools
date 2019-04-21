#!/usr/bin/env python

from rna_tools.BlastPDB import BlastPDB
from rna_tools.RfamSearch import RNASequence, RfamSearch
from rna_tools.rna_tools_lib import RNAStructure


def test_blastpdb():
    p = BlastPDB('GGGUCAGGCCGGCGAAAGUCGCCACAGUUUGGGGAAAGCUGUGCAGCCUGUAACCCCCCCACGAAAGUGGG')
    p.search()
    assert p.result.startswith('<HTML') is True


def test_rfamsearch():
    seq = RNASequence("GGCGCGGCACCGUCCGCGGAACAAACGG")
    rs = RfamSearch()
    hit = rs.cmscan(seq)
    print(hit)


def test_rnastructre():
    pass
