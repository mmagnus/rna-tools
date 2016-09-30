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

from utils import *

def parse_args():
    """setup program options parsing"""
    parser = OptionParser()
    parser.add_option("--data-dir", dest="data_dir",
                  help="data dir", metavar="DIR", default=DATA_DIR)
    (options, args)  = parser.parse_args()

    return (parser, options, args)

def show_doublet_info(doublet_id, options):
    pdb_id,d1,d2 = doublet_id.split(":")
    pdb_id = pdb_id.upper()

    full_fn = PDBObject.pdb_fn(pdb_id,"res_map",data_dir=options.data_dir)
    print "info from %s" % os.path.basename(full_fn)
    if not os.path.isfile(full_fn):
        print " -- missing file --"
    else:
        d_dict = load_json(full_fn)
        (c1,num1,c2,num2) = (d1[0],d1[1:],d2[0],d2[1:])
        org_info1 = d_dict.get(c1,{}).get(num1,{})
        org_num1 = org_info1.get('resseq','???') + org_info1.get('icode','')
        org_info2 = d_dict.get(c2,{}).get(num2,{})
        org_num2 = org_info2.get('resseq','???') + org_info2.get('icode','')
        print " original number: %s:%s%s:%s%s" % (pdb_id,c1,org_num1,c2,org_num2)

    print "info in graph files"
    for fn in ('close_doublets','contacts_RV','contacts_MC','contacts_FR','contacts_MO','contacts_CL'):
        full_fn = PDBObject.pdb_fn(pdb_id,fn,data_dir=options.data_dir)
        print "info from %s" % fn
        if not os.path.isfile(full_fn):
            print " -- missing file --"
            continue
        g = load_graph(full_fn)
        edges = g.get_edge_data(d1,d2)
        if edges is not None:
            for edge in edges.values():
                print edge
        else:
            print " -- no data --"
    print "groups:"
    i = 0
    groups = load_json(PDBObject.pdb_fn(pdb_id,"groups",data_dir=options.data_dir))
    for k in sorted(groups.keys()):
        if doublet_id in groups[k]:
            print " in %s" % k
            i += 1
    if i==0:
        print " -- no data --"

def main():
    (parser, options, args) = parse_args()
    if len(args)!=1:
        print "give doublet id"
        parser.print_help()
        exit(1)
    show_doublet_info(args[0], options)

if __name__ == '__main__':
    main()

