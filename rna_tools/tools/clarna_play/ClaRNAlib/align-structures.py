#!/usr/bin/env python
import sys
import os
import re
import copy
from optparse import OptionParser
from normalize_pdb import normalize_pdb
from itertools import combinations
from scipy.spatial import KDTree

from utils import *

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""aligns input PDB to the reference structure""")
    parser.add_option("-o", "--output", dest="output",
                  help="save result to file", metavar="FILE",
                  default="contacts-db")
    parser.add_option("-i", "--input", dest="input",
                  help="read structure to be aligned", metavar="FILE")
    parser.add_option("-r", "--reference", dest="ref",
                  help="read reference structure", metavar="FILE")
    parser.add_option("--each-model", dest="each_model", action="store_true",
                  help="align each model separately", 
                  default=False)
    (options, args)  = parser.parse_args()
    return (parser, options, args)

def compute_transformation(c_ref,c2):
    sup = SVDSuperimposer()
    sup.set(c_ref, c2)
    sup.run()
    rms = sup.get_rms()    
    (rot,tran) = sup.get_rotran()
    return (rms,rot,tran)

def get_points(s):
    atoms_list = ["C1'","C4","C6","C2"]
    res = []
    if hasattr(s, 'get_chains'):
        chains = s.get_chains()
    else:
        chains = [c for c in s]
    for c_num,c in enumerate(chains):
        for r in c:
            rr = simplify_residue(r)
            if not rr.has_key("C1'"):
                continue
            for a in atoms_list:
                if rr.has_key(a):
                    res.append(rr[a])
                else:
                    return None
                    # res.append(None)
    if len(res)!=2*len(atoms_list):
        return None
    return array(res,'f')
      
def transform_structure(s,transform):
    (rot, tran) = transform
    for a in s.get_atoms():
        a.transform(rot, tran)    

def align_structures(ref_fn, inp_fn, options):
    parser = PDB.PDBParser()
    ref_model = parser.get_structure('ref',ref_fn)[0] # first model
    inp_struct = parser.get_structure('inp', inp_fn)

    ref_points = get_points(ref_model)
    inp_points = get_points(inp_struct[0])
    (rms,rot,tran) = compute_transformation(ref_points, inp_points)
    if options.each_model:
        for m in inp_struct:
            inp_points = get_points(m)
            (rms,rot,tran) = compute_transformation(ref_points, inp_points)
            print "rms=%.5f" % rms
            transform_structure(m, (rot,tran))
    else:
        print "rms=%.5f" % rms
        transform_structure(inp_struct, (rot,tran))

    fn = options.output
    print "saving: %s" % fn
    io = PDB.PDBIO()
    io.set_structure(inp_struct)
    io.save(fn)



def main():
    (parser, options, args) = parse_args()
    align_structures(options.ref, options.input, options)


if __name__ == '__main__':
    main()

