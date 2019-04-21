#!/usr/bin/env python

"""
Example::

    ./rna_helix_vis.py -p test_data/rp14_farna_eloop.out.1.pdb -s test_data/rp14.ss
    (((((.(((((....)))))..[.....((((((..........))))))....)))))].

    /// copy below to pymol ///
    from pymol.cgo import *
    from pymol import cmd
    color red, resi 1+2+3+4+5+55+56+57+58+59
    color red, resi 7+8+9+10+11+16+17+18+19+20
    color red, resi 29+30+31+32+33+34+45+46+47+48+49+50
    obj = [ CYLINDER, -4.89364,-21.9653,-24.0907,10.3337,-22.0812,-28.8867, 9, 1, 0, 0, 1,0,0, CYLINDER, 9.52677,-27.707,-35.113,4.12388,-23.5705,-48.3105, 9, 1, 0, 0, 1,0,0, CYLINDER, 3.48367,-49.4752,-18.6601,-12.6687,-57.1961,-17.8483, 9, 1, 0, 0, 1,0,0,]
    cmd.load_cgo(obj,'cgo01')
    /// copy above to pymol ///

.. image:: ../../rna_tools/tools/rna_helix_vis/doc/rna_helix_vis.png

.. warning:: works only for chain A @todo

required: forgi (https://github.com/pkerpedjiev/forgi)

`pip install forgi`
"""
import argparse
import sys
import six

from rna_tools.tools.rna_bp.base_pair import *

if six.PY2:
    import forgi.graph.bulge_graph as fgb
    bg = fgb.BulgeGraph()
else:
    print('Run under Python 2')


class Helix:
    def __init__(self, id, fn, chain, start, end):
        self.id = id
        self.start = start  # (1,2)
        self.end = end
        self.fn = fn
        self.chain = chain
        self.load_structure()

    def __repr__(self):
        return 'Helix ' + str(self.id) + ' ' + str(self.start) + '-' + str(self.end)

    def get_ids(self):
        t = ''
        for x in range(self.start[0], self.start[1] + 1):
            t += str(x) + '+'
        for x in range(self.end[0], self.end[1] + 1):
            t += str(x) + '+'
        return 'color red, resi ' + t[:-1]

    def load_structure(self):
        parser_pdb = PDB.PDBParser()
        struct = parser_pdb.get_structure('', self.fn)  # )
        model = struct[0]
        self.chain = model[self.chain]
        self.struct = struct

    def get_helix_coord(self):
        chain = self.chain
        struct = self.struct
        a = chain[self.start[0]]
        b = chain[self.end[1]]

        bp = BasePair(a, b, struct)

        a = chain[self.start[1]]
        b = chain[self.end[0]]

        bp2 = BasePair(a, b, struct)

        # print "cmd.load_cgo( [9.0, " + ','.join([str(x) for x in bp.calc_coord()]) + ',' +  ','.join([str(x) for x in bp2.calc_coord()]) + ", 5, 1, 0, 0, 1,0,0], 'cylinderx' )"
        return ','.join([str(x) for x in bp.calc_coord()]) + ',' + ','.join([str(x) for x in bp2.calc_coord()]) + ", 9, 1, 0, 0, 1,0,0"


def show_cyliders(ss):

    bg = fgb.BulgeGraph()
    bg.from_dotbracket(ss)

    helices = []

    print("""from pymol.cgo import *
from pymol import cmd""")

    ht = 'obj = [ CYLINDER, '
    for i in bg.to_bg_string().split('\n'):
        # print i
        if i.startswith('define s'):
            x = i.replace('define s', '')
            # print x
            x = [int(x) for x in x.split()]
            h = Helix(x[0], fn, chain, (x[1], x[2]), (x[3], x[4]))
            helices.append(h)
            print((h.get_ids()))
            ht += h.get_helix_coord() + ', CYLINDER, '

    print((ht[:-11] + ']'))
    print("cmd.load_cgo(obj,'cgo01') ")


def get_parser():
    parser = argparse.ArgumentParser()  # usage="%prog [<options>] <pdb files (test_data/*)>")
    parser.add_argument('-p', "--pdb_fn",
                        dest="pdb_fn",
                        default='',
                        help="pdb file")

    parser.add_argument('-s', "--ss_fn",
                        dest='ss',
                        default='',
                        help="secondary structure file, only with ss (one line, one chain)")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    chain = 'A'
    fn = args.pdb_fn
    ss = open(args.ss).read().strip()

    print(ss)
    print('/// copy below to pymol ///')
    show_cyliders(ss)
    print('/// copy above to pymol ///')
