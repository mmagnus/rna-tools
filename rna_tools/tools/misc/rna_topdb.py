#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
                                   PREDS                                     REAL
^^--> 0 OP2 res C group 0 compare  tensor([29.6149, 36.8916, 42.6567])  in:  tensor([29.6170, 36.8990, 42.6580]) RMSD  0.007772309705796946

input:

ATOM     20  OP2   C A  18      29.617  36.899  42.658  1.00 75.10           O1-

requried:

   pip install rna-tools

"""
import argparse
import re


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("input", help="output file", default="") # nargs='+')
    parser.add_argument("output", help="pdb file", default="")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    from rna_tools.rna_tools_lib import RNAStructure
    
    rna = RNAStructure()
    for l in open(args.input):
            t = "\^\^\-\-\> (?P<resid>\d+) (?P<atomname>[\w']+) res (?P<res>\w+).+tensor\(\[(?P<xyz>.+)\]\)"
            r = re.search(t, l)
            if r:
                print(l.strip())
                
                resid = int(r.group('resid')) + 1
                atomname = r.group('atomname')
                res = r.group('res')
                x,y,z = [float(x) for x in r.group('xyz').split(',')]
                
                l = rna.get_empty_line()
                
                l = rna.set_res_code(l, res)
                l = rna.set_atom_code(l, atomname)
                l = rna.set_res_index(l, resid)
                l = rna.set_atom_coords(l, x, y, z)
                
                rna.add_line(l)
                #print(resid, atomname, res, x, y, z)
        
    rna.write(args.output)
