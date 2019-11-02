#!/usr/bin/env python

from __future__ import print_function
import os
import sys

from rna_tools.BlastPDB import BlastPDB
from rna_tools.RfamSearch import RNASequence, RfamSearch
from rna_tools.rna_tools_lib import RNAStructure
from rna_tools.rna_tools_lib import edit_pdb, add_header, get_version, \
                          collapsed_view, fetch, fetch_ba, replace_chain, RNAStructure, \
                          select_pdb_fragment, sort_strings
from rna_tools.tools.rna_x3dna.rna_x3dna import x3DNA


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

if __name__ == '__main__':
    # get version
    version = os.path.basename(os.path.dirname(os.path.abspath(__file__))), get_version(__file__)
    version = version[1].strip()

    print('- Python:', sys.version.replace('\n', ''))

    print('- rna-tools:', version)

    print('- RNA_TOOLS_PATH set to ', os.environ['RNA_TOOLS_PATH'])

    print('- See full list of tools <https://github.com/mmagnus/rna-tools/blob/master/rna-tools-index.csv')
    
    print('Seems OK')




    
