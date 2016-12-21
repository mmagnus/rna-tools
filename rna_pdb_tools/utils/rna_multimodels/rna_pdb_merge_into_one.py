#!/usr/bin/env python

import sys
from rna_pdb_tools.pdb_parser_lib import StrucFile

if __name__ == '__main__':
    files = sys.argv[1:]
    if not files:
        print 'rna-pdb-merge-into-one.py test_in/*.pdb'
        sys.exit(1)

    c = 1
    for f in files:
        s = StrucFile(f)
        print ("MODEL        " + str(c))
        print s.get_text(add_end=False)
        print ('ENDMDL')
        c += 1
    print ('END')
