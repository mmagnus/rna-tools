#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert a secondary structure with a pk to the SimRNA format::

       rna_convert_pseudoknot_formats git:(master) ✗ python rna_ss_pk_to_simrna.py test_data/ss_with_pk.ss
       ((((([[[[[[)))))........(.((....(]]]]]].)..(((. .)))...)).)

       (((((......)))))........(.((....(.......)..(((. .)))...)).)
       .....((((((......................))))))........ ...........
"""

import argparse

SECOND_PK_CHAR = ['{', '}']

def is_pk(ss):
    if ss.find('[') > -1:
        return True
    if ss.find(']') > -1:
        return True

def get_multiple_lines(ss):
    ss_n = ''
    ss_pk = ''
    for i in ss.strip():
        if i == '[':
            j = '('
            i = '.'
        elif i == ']':
            j = ')'
            i = '.'
        elif i == ' ': # not space
            j = ' '
        else:
            j = '.'

        ss_n += i
        ss_pk += j

    return ss_n + '\n' + ss_pk

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='file')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    if args.file:
        ss = open(args.file).readline()
        print(ss)
        if is_pk(ss):
            print((get_multiple_lines(ss)))
        else:
            print(ss)
