#!/usr/bin/env python
import sys
import os
import re
import copy
from optparse import OptionParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

from utils import *
from distances import NORMALIZED_BASE

def parse_args():
    """setup program options parsing"""
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="output",
                  help="save result to file", metavar="FILE",
                  default="contacts-db")
    parser.add_option("--data-dir", dest="data_dir",
                  help="directory with data", metavar="DIR", default="gc-data")
    parser.add_option("--doublet-id", dest="doublet_id",
                  help="read structure to be aligned", metavar="PDBID:R1:R2")
    parser.add_option("--chain-a", dest="chain_a",
                  help="replace chain A with given nucleotide", metavar="PDBID:R1|RESNAME")
    parser.add_option("--chain-b", dest="chain_b",
                  help="replace chain B with given nucleotide", metavar="PDBID:R1|RESNAME")
    (options, args)  = parser.parse_args()
    return (parser, options, args)

def compute_transformation(c_ref,c):
    sup = SVDSuperimposer()
    sup.set(c_ref, c)
    sup.run()
    rms = sup.get_rms()
    (rot,tran) = sup.get_rotran()
    return (rms,rot,tran)

def morph_residue(org_points, org_name, new_points, new_name):
    # assert sorted([org_name, new_name]) in (['C','C'],['C','U'],['U','U'],['A','A'],['A','G'],['G','G'])
    if org_name in ['C','U']:
        A1 = ['C6','C4','C2']
    else:
        A1 = ['C4','C6','C2']
    if new_name in ['C','U']:
        A2 = ['C6','C4','C2']
    else:
        A2 = ['C4','C6','C2']
    # A1,A2 = ['C2','C4','C6'],['C2','C4','C6']
    org_arr = array([org_points[a] for a in A1],'f')
    new_arr = array([new_points[a] for a in A2],'f')
    (rms,rot,tran) = compute_transformation(org_arr, new_arr)
    result = {}
    for key,p in new_points.items():
        result[key] = np.dot(np.array(p,'f'),rot)+tran
    return result

def morph(options):
    dd = DoubletsDict(options.data_dir, reduced_atoms=['*'])

    (p1,p2) = dd.get_normalized(options.doublet_id)
    n_type = dd.get_n_type(options.doublet_id)
    
    org_n1,org_n2 = n_type[0],n_type[1]
    new_n1,new_n2 = n_type[0],n_type[1]
    
    
    if options.chain_a:
        if re.match("^[a-zA-Z0-9]*:[a-zA-Z0-9]*$",options.chain_a):
            new_n1 = dd.get_residue_name(options.chain_a)
            new_p1 = dd.get_residue(options.chain_a)
        else:
            new_n1 = options.chain_a.upper()
            new_p1 = NORMALIZED_BASE[new_n1]
        print "morphing chain-A with residue: %s (%s->%s)" % (options.chain_a, org_n1, new_n1)
        p1 = morph_residue(p1, n_type[0], new_p1, new_n1)

    if options.chain_b:
        if re.match("^[a-zA-Z0-9]*:[a-zA-Z0-9]*$",options.chain_b):
            new_n2 = dd.get_residue_name(options.chain_b)
            new_p2 = dd.get_residue(options.chain_b)
        else:
            new_n2 = options.chain_b.upper()
            new_p2 = NORMALIZED_BASE[new_n2]
        print "morphing chain-B with residue: %s (%s->%s)" % (options.chain_b, org_n2, new_n2)
        p2 = morph_residue(p2, n_type[1], new_p2, new_n2)

    new_n_type = new_n1+new_n2

    gp = __import__('gen-pdb')
    s = PDB.Structure.Structure("points")
    for ignored in ["NEXT:O3'"]:
        if p1.has_key(ignored):
            del p1[ignored]
        if p2.has_key(ignored):
            del p2[ignored]
    s.add(gp.points2model(1,(p1,p2),new_n_type))
    save_pdb(options.output,s)


def main():
    (_parser, options, _args) = parse_args()
    morph(options)

if __name__ == '__main__':
    main()

