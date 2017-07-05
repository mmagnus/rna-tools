#!/usr/bin/env python

try:
    import urllib2
except:
    print('BlastPDB requires urllib2')

class BlastPDB:
    """BlastPDB - run Blast online on the PDB database.

    This can be used in Jupiter based RNA notebooks, e.g. 
    https://github.com/mmagnus/rna-pdb-tools/blob/master/rp18.ipynb
    
    Usage::

       >>> p = BlastPDB('GGGUCAGGCCGGCGAAAGUCGCCACAGUUUGGGGAAAGCUGUGCAGCCUGUAACCCCCCCACGAAAGUGGG')
       >>> p.search()
       >>> p.result  #doctest: +ELLIPSIS
       '<HTML>\\n<TITLE>BLAST Search Results</TITLE>...

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

    import doctest
    doctest.testmod()
