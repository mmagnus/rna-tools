#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Secondary structure analysis"""

import os
import tempfile
import shutil
import subprocess

from rna_tools.rna_tools_config import VARNA_JAR_NAME, VARNA_PATH


class ExceptionOpenPairsProblem(Exception):
    pass


def draw_ss(title, seq, ss, img_out, resolution=4, verbose=False):
    """Draw Secondary Structure using VARNA (you need correct configuration for this).

    If everything is OK, return None, if an error (=exception) return stderr.

    Usage::

       >>> seq = 'GGAAACC'
       >>> ss =  '((...))'
       >>> img_out = 'output/demo.png'
       >>> draw_ss('rna', seq, ss, img_out)
       >>> print('Made %s' % img_out)
       Made output/demo.png

    .. image:: ../../rna_tools/output/demo.png
       :scale: 25 %

    Can be used with http://geekbook.readthedocs.io/en/latest/rna.html"""
    curr = os.getcwd()
    os.chdir(VARNA_PATH)  # VARNAv3-93-src')
    if verbose:
        print(VARNA_PATH)
    t = tempfile.NamedTemporaryFile(delete=False)
    t.name += '.png'

    cmd = 'java -cp ' + VARNA_JAR_NAME + ' fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN ' + seq + \
        " -structureDBN '" + ss + "' -o " + t.name + " -title '" + \
        title + "' -resolution '" + str(resolution) + "'"
    if verbose:
        print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out = p.stderr.read().decode().strip()
    os.chdir(curr)
    if out.find('Exception') > -1:
        return out
    else:
        if verbose:
            print(t.name)
        shutil.move(t.name, img_out)


def parse_vienna_to_pairs(ss, remove_gaps_in_ss=False):
    """Parse Vienna (dot-bracket notation) to get pairs.

    Args:

       ss (str): secondary stucture in Vienna (dot-bracket notation) notation
       remove_gaps_in_ss (bool): remove - from ss or not, design for DCA (tpp case
                                 ``ss = "(((((((((.((((.(((.....))))))......------)....."``
                                 works with pk of the first level, ``[[]]``

    Returns:

        list of two lists: (pairs, pairs_pk)

    Examples::

        >>> parse_vienna_to_pairs('((..))')
        ([[1, 6], [2, 5]], [])

        >>> parse_vienna_to_pairs('(([[))]]')
        ([[1, 6], [2, 5]], [[3, 8], [4, 7]])

        >>> parse_vienna_to_pairs('((--))')
        ([[1, 6], [2, 5]], [])

        >>> parse_vienna_to_pairs('((--))', remove_gaps_in_ss=True)
        ([[1, 4], [2, 3]], [])

        >>> parse_vienna_to_pairs('((((......')
        Traceback (most recent call last):
          File "/usr/lib/python2.7/doctest.py", line 1315, in __run
            compileflags, 1) in test.globs
          File "<doctest __main__.parse_vienna_to_pairs[4]>", line 1, in <module>
            parse_vienna_to_pairs('((((......')
          File "./SecondaryStructure.py", line 106, in parse_vienna_to_pairs
            raise ExceptionOpenPairsProblem('Too many open pairs (()) in structure')
        ExceptionOpenPairsProblem: Too many open pairs (()) in structure

    """
    if remove_gaps_in_ss:
        ss = ss.replace('-', '')
    stack = []
    pairs = []
    pairs_pk = []
    stack_pk = []
    for c, s in enumerate(ss):
        if s == '(':
            stack.append(c + 1)
        if s == ')':
            pairs.append([stack.pop(), c + 1])
        if s == '[':
            stack_pk.append(c + 1)
        if s == ']':
            pairs_pk.append([stack_pk.pop(), c + 1])

    if stack:
        raise ExceptionOpenPairsProblem('Too many open pairs (()) in structure')
    if stack_pk:
        raise ExceptionOpenPairsProblem('Too many open pairs [[]] in structure')

    pairs.sort()
    pairs_pk.sort()
    return(pairs, pairs_pk)


# main
if __name__ == '__main__':
    import doctest
    doctest.testmod()
