#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Set of functions for working with single lines from PDB files
"""

def get_atom_num(line):
    """Extract atom number from a line of PDB file

    Arguments:
      * line = ATOM line from a PDB file

    Output:
      * atom number as an integer
    """
    return int(''.join([x for x in line[6:10] if x.isdigit()]))


def get_res_num(line):
    """Extract residue number from a line of PDB file

    Arguments:
      * line = ATOM line from a PDB file

    Output:
      * residue number as an integer
    """
    return int(''.join([x for x in line[22:27] if x.isdigit()]))


def get_res_code(line):
    """Get residue code from a line of a PDB file
    """
    if not line.startswith('ATOM'):
        return None
    return line[17:20]


def get_atom_code(line):
    """Get atom code from a line of a PDB file
    """
    if not line.startswith('ATOM'):
        return None
    return line[13:16].replace(' ', '')

def get_atom_coords(line):
    """Get atom coordinates from a line of a PDB file
    """
    if not line.startswith('ATOM'):
        return None
    return tuple(map(float, line[31:54].split()))


def set_line_bfactor(line, bfactor):
    if not line.startswith('ATOM'):
        return None
    return line[:60] + (" %5.2f" % bfactor) + line[66:]


def set_atom_code(line, code):
    return line[:13] + code + ' ' * (3 - len(code)) + line[16:]


def set_res_code(line, code):
    return line[:18] + code.rjust(3) + line[21:]
