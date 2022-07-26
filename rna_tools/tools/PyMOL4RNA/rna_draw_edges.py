#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import argparse
from Bio import PDB

    
def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--name', default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    parser = PDB.PDBParser()
    struct = parser.get_structure('', args.file)
    model = struct[0]
    print("cmd.delete('%s_edges')" % args.name)
    print("cmd.do('radius=0.1')")
    for chain in model:
        for r in chain:
            chain = r.get_parent().id
            resi = r.get_id()[1]

            # this is sooooooo ugly #
            if r.get_resname().strip().startswith('U'):
                x1,y1,z1 = [str(x) for x in r['O4'].get_coord()]
                x2,y2,z2 = [str(x) for x in r['O2'].get_coord()]
                r1,g1,b1,r2,g2,b2 = [0.0,0.5,1.0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x2 + ',' + y2 +',' + z2 +', radius, %i,%i,%i,%i,%i,%i], "WatsonCrick_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))
                x3,y3,z3 = [str(x) for x in r["C5"].get_coord()]
                r1,g1,b1,r2,g2,b2 = [1,0,0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x3 + ',' + y3 +',' + z3 +', radius,  %i,%i,%i,%i,%i,%i], "Hoogesteen_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))
                atom = r["O2'"]
                x4,y4,z4 = [str(x) for x in atom.get_coord()]
                r1,g1,b1,r2,g2,b2 = 1,1,0,1,1,0
                print('cmd.load_cgo( [ 9.0, ' + x2 + ','+ y2 + ','+ z2 +',' + x4 + ',' + y4 +',' + z4 +', radius,  %i,%i,%i,%i,%i,%i], "Sugar_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))

            if r.get_resname().strip().startswith('A'):
                x1,y1,z1 = [str(x) for x in r['N6'].get_coord()]
                x2,y2,z2 = [str(x) for x in r['C2'].get_coord()]
                r1,g1,b1,r2,g2,b2 = [0.0,0.0,1.0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x2 + ',' + y2 +',' + z2 +', radius, %i,%i,%i,%i,%i,%i], "WatsonCrick_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))
                atom = r["N7"]
                x3,y3,z3 = [str(x) for x in atom.get_coord()]
                r1,g1,b1,r2,g2,b2 = [1,0,0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x3 + ',' + y3 +',' + z3 +', radius,  %i,%i,%i,%i,%i,%i], "Hoogesteen_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))
                atom = r["O2'"]
                x4,y4,z4 = [str(x) for x in atom.get_coord()]
                r1,g1,b1,r2,g2,b2 = 1,1,0,1,1,0
                print('cmd.load_cgo( [ 9.0, ' + x2 + ','+ y2 + ','+ z2 +',' + x4 + ',' + y4 +',' + z4 +', radius,  %i,%i,%i,%i,%i,%i], "Sugar_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))

            if r.get_resname().strip().startswith('G'):
                x1,y1,z1 = [str(x) for x in r['O6'].get_coord()]
                x2,y2,z2 = [str(x) for x in r['N2'].get_coord()]
                r1,g1,b1,r2,g2,b2 = [0.0,0.0,1.0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x2 + ',' + y2 +',' + z2 +', radius, %i,%i,%i,%i,%i,%i], "WatsonCrick_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))

                x3,y3,z3 = [str(x) for x in r["N7"].get_coord()]
                r1,g1,b1,r2,g2,b2 = [1,0,0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x3 + ',' + y3 +',' + z3 +', radius,  %i,%i,%i,%i,%i,%i], "Hoogesteen_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))
                atom = r["O2'"]
                x4,y4,z4 = [str(x) for x in atom.get_coord()]
                r1,g1,b1,r2,g2,b2 = 1,1,0,1,1,0
                print('cmd.load_cgo( [ 9.0, ' + x2 + ','+ y2 + ','+ z2 +',' + x4 + ',' + y4 +',' + z4 +', radius,  %i,%i,%i,%i,%i,%i], "Sugar_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))

            if r.get_resname().strip().startswith('C'):
                x1,y1,z1 = [str(x) for x in r['N4'].get_coord()]
                x2,y2,z2 = [str(x) for x in r['O2'].get_coord()]
                r1,g1,b1,r2,g2,b2 = [0.0,0.0,1.0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x2 + ',' + y2 +',' + z2 +', radius, %i,%i,%i,%i,%i,%i], "WatsonCrick_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))

                x3,y3,z3 = [str(x) for x in r["C6"].get_coord()]
                r1,g1,b1,r2,g2,b2 = [1,0,0] * 2
                print('cmd.load_cgo( [ 9.0, ' + x1 + ','+ y1 + ','+ z1 +',' + x3 + ',' + y3 +',' + z3 +', radius,  %i,%i,%i,%i,%i,%i], "Hoogesteen_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))
                atom = r["O2'"]
                x4,y4,z4 = [str(x) for x in atom.get_coord()]
                r1,g1,b1,r2,g2,b2 = 1,1,0,1,1,0
                print('cmd.load_cgo( [ 9.0, ' + x2 + ','+ y2 + ','+ z2 +',' + x4 + ',' + y4 +',' + z4 +', radius,  %i,%i,%i,%i,%i,%i], "Sugar_%s-%i" )' % (r1,g1,b1,r2,g2,b2, chain, resi))

    print("cmd.group('%s_edges','WatsonCrick*')" % args.name)
    print("cmd.group('%s_edges','Hoogesteen_*')" % args.name)
    print("cmd.group('%s_edges','Sugar_*')" % args.name)
