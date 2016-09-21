#!/usr/bin/env python

from collections import OrderedDict
import re
import string
import sys

def select_pdb_fragment_pymol_style(txt):
    """Take txt such as A/10-15/P and parse into::
    
    A/57/O2' -> ['A', ['57'], "O2'"]

    .. warning:: e.g. for A:1-31, resi 31 is not included, works like in Python."""
    v = 0
    selection = OrderedDict()
    # 
    if txt.find('|') < -1:
        txt = txt.replace(' ','').split('|')
    else:
        txt = [txt]

    for t in txt:
        l = t.split('/')
        if v:print l
        l[1] = l[1].split('+')

        if l[0] in string.ascii_letters:
            chain_curr = l[0]

        for i in l[1]:
            if i.find('-') > -1:
                start, ends = i.split('-')
                if start > ends:
                    print >>sys.stderr, 'Error: range start > end ' + i
                    return False
                index = range(int(start), int(ends))#+1)
            else:
                index=[int(i)]

            index_and_atoms = [index, l[2].split('+')]
            if v: print index_and_atoms

            if selection.has_key(chain_curr):
                selection[chain_curr] += [index_and_atoms]
            else:
                selection[chain_curr] = [index_and_atoms]

    return selection

def is_in_selection(selection, curr_chain_id, curr_resi, curr_atom_name):
    if selection.has_key(curr_chain_id):
            for sele_range in selection[curr_chain_id]:
                if curr_resi in sele_range[0]:
                        if curr_atom_name in sele_range[1]:
                            return True
    return False                  

#main
if __name__ == '__main__':
    selection = select_pdb_fragment_pymol_style('E/1-15+18/P | A/1-3/P')

    curr_chain_id = 'E'
    curr_resi = 1
    curr_atom_name = 'P'

    print selection
    if selection:
        print is_in_selection(selection, curr_chain_id, curr_resi, curr_atom_name)

    print is_in_selection(selection, 'X', curr_resi, curr_atom_name)        
    print is_in_selection(selection, 'E', curr_resi, "C'")        
    print is_in_selection(selection, 'A', curr_resi, "P'")        
