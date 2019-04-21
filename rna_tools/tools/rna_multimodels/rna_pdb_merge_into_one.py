#!/usr/bin/env python

import sys
from rna_tools.rna_tools_lib import RNAStructure


if __name__ == '__main__':
    files = sys.argv[1:]
    if not files:
        print('rna_pdb_merge_into_one.py test_in/*.pdb')
        sys.exit(1)

    c = 1
    for f in files:
        s = RNAStructure(f)
        print("MODEL        " + str(c))
        print(s.get_text(add_end=False))
        print('ENDMDL')
        c += 1
    print('END')
