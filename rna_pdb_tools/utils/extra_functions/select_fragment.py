#!/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
import re
import string
import sys

def select_pdb_fragment(txt, separator="-", splitting='[:\+]', verbose=False):
    """Take txt such as ``A:1-31+B:1-11`` and parse into::

          OrderedDict([('A', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
          15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]),
          ('B', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])])          
    
      .. warning:: e.g. for A:1-31, resi 31 is included"""
    txt = txt.replace(' ','')
    if verbose: print(txt)
    #l = re.split, txt)
    l = re.split(splitting, txt)
    if verbose: print(l)

    selection = OrderedDict()
    for i in l: # ['A', '1-10', '15', '25-30', 'B', '1-10']
        if i in string.ascii_letters:
            if verbose: print('chain', i)
            chain_curr = i
            continue

        if i.find(separator) > -1:
            start, ends = i.split(separator)
            start = int(start)
            ends = int(ends)
            if start > ends:
                print('Error: range start > end ' + i, file=sys.stderr)
                sys.exit(1)
            index = list(range(int(start), int(ends)+1)) # without +1 python like, with +1 people-like
        else:
            index=[int(i)]
        if chain_curr in selection:
            selection[chain_curr] += index
        else:
            selection[chain_curr] = index
    return selection


def select_pdb_fragment_pymol_style(txt):
    """Take txt such as A/10-15/P and parse into::
    
      A/57/O2' -> ['A', ['57'], "O2'"]

    If you want to combine a few subselections, please use ``,``::
 
       --model_ignore_selection "A/57/O2',A/58/O2'"
    
    .. warning:: e.g. for A:1-31, resi 31 is included"""
    v = 0
    selection = OrderedDict()
    if txt.find(',') > -1:
        txt = txt.replace(' ','').split(',')
    else:
        txt = [txt]
    for t in txt:
        l = t.split('/')
        if v:print(l)
        l[1] = l[1].split('+')

        if l[0] in string.ascii_letters:
            chain_curr = l[0]

        for i in l[1]:
            if i.find('-') > -1:
                start, ends = i.split('-')
                if start > ends:
                    print('Error: range start > end ' + i, file=sys.stderr)
                    return False
                index = list(range(int(start), int(ends)))#+1)
            else:
                index=[int(i)]

            index_and_atoms = [index, l[2].split('+')]
            if v: print(index_and_atoms)

            if chain_curr in selection:
                selection[chain_curr] += [index_and_atoms]
            else:
                selection[chain_curr] = [index_and_atoms]

    return selection

def is_in_selection(selection, curr_chain_id, curr_resi, curr_atom_name):
    if curr_chain_id in selection:
            for sele_range in selection[curr_chain_id]:
                if curr_resi in sele_range[0]:
                        if curr_atom_name in sele_range[1]:
                            return True
    return False                  

#main
if __name__ == '__main__':
    selection = select_pdb_fragment_pymol_style('E/1-15+18/P, A/1-3/P')

    curr_chain_id = 'E'
    curr_resi = 1
    curr_atom_name = 'P'

    print(selection)
    if selection:
        print(is_in_selection(selection, curr_chain_id, curr_resi, curr_atom_name))

    print(is_in_selection(selection, 'X', curr_resi, curr_atom_name))        
    print(is_in_selection(selection, 'E', curr_resi, "C'"))        
    print(is_in_selection(selection, 'A', curr_resi, "P'"))        

    print(select_pdb_fragment_pymol_style('A/48/OP2,B/48/OP2'))

    print(select_pdb_fragment('A:1-31+B:1-11'))
