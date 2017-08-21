#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This beauty here will go to sali notation and convert it to dotbracket notation.
The file name should be xxxx.sali

Author: Catarina Almeida"""

import re
import string
import sys
import argparse
import os

dic_ss_elements = {'>': ')', '<': '(', '[': '(', ']': ')', '~': '-', '^': '.'}


def repl(m):
    """This function will replace the length of a given string by the correspondent number of dashes.
    The expression ``qwerty`` will be replaced by ``-----``.
    """
    return '-' * len(m.group())


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help="file in the Sali format")
    return parser


def convert_sali2dotbracket(fn):
    """The function needs a filename in the Sali format.
    This function will get the secondary structure of the sequence,
    then get its identifier and then the sequence itself.

    **To get the ss**

    The line with the secondary structure is a list and will look like this::

       ['', '', '', '', '', '', '', '', '', '', '--...<<<[[...]]..>>>>', '', '', '\\n']

    In this case, the ss is in the 11th position. But in some files it may be in the 12th, 13th, 10th, etc..

    If the longest element from the list is extracted, then this problem is overcomed.

    The ss will some times have patterns of repeated gaps, which will come in the form of:

    a. x

    b. xnt

    c. ( x )

    With x being any number, from 1 to 1000.
    These must be converted to the correspondent number of gaps (-) in the converted ss.
    This conversion is done by:

    1 - Identifying the pattern with regex

    2 - Replacing it with repl function.

    As such, the following expressions will replace the previously mentioned patterns:

    a. ``re.sub(r'\d*\d', repl, temp)``

    b. ``re.sub(r'\d*\dnt', repl, temp)``

    c. ``re.sub(r'(?P<smthBeautiful>\(\d+\))', repl, temp)``

    **To get the sequence**

    The sequence, much like the ss, can sometimes be in a different position in the list. Like in the ss, the longest element will be selected.
    Also, like in the ss, patterns for repeated gaps appear. So these must also be removed."""
    for line in open(fn):
        if line.startswith(' ') and not re.search('[a-zA-Z]', line):
            # see pydoc # ss #
            lineSplit = line.split(' ')
            ssbefore = max(lineSplit, key=len)
            ssafter = ''
            for x in ssbefore:
                if x in dic_ss_elements:
                    ssafter += dic_ss_elements[x]
                else:
                    ssafter += x
            re.sub(r'(?P<smthBeautiful>\(\d+\))', repl, ssafter)  # (951)
            re.sub(r'\d*\dnt', repl, ssafter)  # 1598nt
            print(ssafter)
        elif line.startswith(' ') or line.startswith('#'):
            True
        else:
            # get the name of the file and add it to the identifier
            nho = os.path.basename(fn).split('.')[0]
            print('>' + nho + '_' + line.split(" ")[0])  # get the identifiers in fasta format
            # see pydoc # sequence #
            temp = line.split(' ')
            sequence = max(temp, key=len)  # the longest element of the list will be the sequence
            print(re.sub(r'\d*\dnt|\d*\d|(?P<smthBeautiful>\(\d+\))', repl, sequence, end=' '))


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    convert_sali2dotbracket(args.filename)
