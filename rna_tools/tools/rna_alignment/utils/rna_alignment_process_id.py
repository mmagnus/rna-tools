#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("-d", "--debug",
                        action="store_true", help="be verbose")
    parser.add_argument('--id-width', type=int, default=70)
    parser.add_argument('--sep', default='|')
    parser.add_argument("alignment")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    name_width = args.id_width

    for l in open(args.alignment):
        if l:
            if args.debug:
                print(l)
            if not l.startswith('#') and not l.startswith('//'):
                if args.debug: print(l.split())
                id, seq = l.split()
                id = id.replace('[', '|').replace(']', '|')
                if args.debug: print(id)
                species, group = ('','')
                if len(id.split('|')) > 2:
                    species, group = id.split('|')[:2]
                    if args.debug: print(group, species)
                    id = group + '|' + species
                else:
                    id = id[0]
                # Leishmania-major-strain-Friedlin[Euglenozoa]FR796420.1|1.0|CUCU-AUG/1-7
                line = id.ljust(name_width) + seq.strip()
                print(line)
            elif '#=GC RF_cons' in l:
                ss = l.replace('#=GC RF_cons', '')
                print('#=GC RF_cons'.ljust(name_width) + ss.strip())
            elif '#=GC SS_cons' in l:
                ss = l.replace('#=GC SS_cons', '')
                print('#=GC SS_cons'.ljust(name_width) + ss.strip())
            else:
                print(l.strip())
