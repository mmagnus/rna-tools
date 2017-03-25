#!/usr/bin/env python

import urllib2

class BlastPDB:
    """BlastPDB

    :param seq: string
    """
    def __init__(self, seq):
        self.seq = seq
        self.result = ''

    def search(self):
        """Search online the seq."""
        p = urllib2.urlopen('http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence=' + self.seq + '&eCutOff=10.0&matrix=BLOSUM62&outputFormat=HTML')
        self.result = p.read()

if __name__ == '__main__':
    p = BlastPDB('GGGUCAGGCCGGCGAAAGUCGCCACAGUUUGGGGAAAGCUGUGCAGCCUGUAACCCCCCCACGAAAGUGGG')
    p.search()
    print(p.result)
