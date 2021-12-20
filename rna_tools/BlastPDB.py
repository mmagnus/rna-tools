#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import urllib3
except:
    print('BlastPDB requires urllib3')


class BlastPDB:
    """BlastPDB - run Blast online on the PDB database.

    This can be used in Jupiter based RNA notebooks, e.g.
    https://github.com/mmagnus/rna-pdb-tools/blob/master/rp18.ipynb

    Warning: getBlastPDB1 has been permanently removed as part of our announced shutdown on December 9th, 2020.
    https://www.rcsb.org/pdb/rest/getBlastPDB1
    
    Usage::

       >>> p = BlastPDB('GGGUCAGGCCGGCGAAAGUCGCCACAGUUUGGGGAAAGCUGUGCAGCCUGUAACCCCCCCACGAAAGUGGG')
       >>> p.search()
       >>> p.result  #doctest: +ELLIPSIS
       u'<HTML>\\n<TITLE>BLAST Search Results</TITLE>...

    :param seq: string
    """

    def __init__(self, seq):
        self.seq = seq
        self.result = ''

    def search(self):
        """Search online the seq."""
        http = urllib3.PoolManager()
        response = http.request('GET', 'http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence=' +
                                self.seq +
                                '&eCutOff=10.0&matrix=BLOSUM62&outputFormat=HTML')
        if response.status == 200:
            self.result = response.data.decode()


# main
if __name__ == '__main__':
    import doctest
    doctest.testmod()
