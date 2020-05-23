#!/usr/bin/env python
import sys
import os
import re
import shutil
import StringIO
import simplejson
import gzip
import tempfile
from optparse import OptionParser
from normalize_pdb import normalize_pdb
from itertools import combinations
from scipy.spatial import KDTree
from progressbar import ProgressBar, Percentage, Bar, ETA

from structure_ciach import StructureCiachCiach
from utils import *

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="create JSON dictionary with residues")
    parser.add_option("-o", "--output-json", dest="output_json",
                  help="save dictionary to JSON file", metavar="DIR")
    parser.add_option("-i", "--input", dest="input",
                  help="parse structures from directory or use given file", metavar="DIR|FILE",
                  default="pdb-simrna")
    parser.add_option("--residues-mapping", dest="residues_mapping",
                  help="parse structures from directory or use given file", metavar="FILE")
    parser.add_option("--with-backbone", dest="with_backbone", action="store_true",
                  help="add backbone fragments of the contacts",
                  default=False)
    parser.add_option("--with-neighbours", dest="with_neighbours", action="store_true",
                  help="add neighbours to the dictionary",
                  default=False)
    parser.add_option("--progress", dest="progress", action="store_true",
                  help="show progress bar",
                  default=False)
    (options, args)  = parser.parse_args()

    return (parser, options, args)


def handle_file(pdb, options):
    if re.match(".*.gz$", pdb):
        f = gzip.open(pdb)
    else:
        f = open(pdb)
    
    if options.residues_mapping:
        r_map = load_json(options.residues_mapping)
    else:
        rm_fn = re.sub(r".pdb(.gz)?$",r".json\1",pdb)
        if os.path.isfile(rm_fn):
            r_map = load_json(rm_fn)
        else:
            print "missing residue mapping dictionary"
            r_map = {}
    
    parser = PDB.PDBParser()
    structure = parser.get_structure("doublet", f)
    if re.match('^.*/zz.*.pdb(.gz)$', pdb):
        req_atoms_list = ('C2','C4','C6')
    else:
        req_atoms_list = ('C2','C4','C6',"C1'",'P',"O3'")
    ciach = StructureCiachCiach(structure,dont_normalize=True,req_atoms_list=req_atoms_list)
    n = len(ciach.good_residues)
    if options.progress:
        widgets = ['Make residues', Percentage(), ' ', Bar(), ' ', ETA()]
        pbar = ProgressBar(widgets=widgets, maxval=n).start()

    res = {}
    for i,id in enumerate(sorted(ciach.good_residues)):
        chain = id[0]
        res_info = r_map.get(chain,{}).get(id[1:],None)
        if res_info is not None:
            org_id = chain+res_info['resseq']+res_info['icode']
        else:
            org_id = None
        
        # TODO: add warning if NEXT:O3' is missing
        structure = ciach.get_single_residue(id,with_backbone=options.with_backbone)
        res[id] = {
            'id': id,
            'org_id': org_id,
            'resname': ciach.get_resname(id),
            'pdb': structure2string(structure),
            'atoms': ciach.get_res_atoms_dict(id),
            'chain': chain,
            'resseq': int(id[1:]),
        }
        if options.with_neighbours:
            res[id]['neighbours'] = ciach.get_neighbours(id)
        if options.progress:
            pbar.update(i+1)
    if re.match(".*.gz$", options.output_json):
        f = gzip.open(options.output_json,"wb")
    else:
        f = open(options.output_json,"w")
    print "saving %d residues to %s" % (len(res.keys()), options.output_json)
    simplejson.dump(res, f)
    if options.progress:
        pbar.finish()



def main():
    (parser, options, args) = parse_args()
    if not options.input:
        print "select input"
        parser.print_help()
        exit(1)
    if not options.output_json:
        print "select output"
        parser.print_help()
        exit(1)
    handle_file(options.input, options)

if __name__ == '__main__':
    main()

