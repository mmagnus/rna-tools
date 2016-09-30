#!/usr/bin/env python
import sys
import os
import re
from utils import *

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""extract contacts from PDB file
""")
    parser.add_option("-o", "--output-json", dest="output_json",
                  help="save doublets to json file", metavar="JSON")
    parser.add_option("-i", "--input-graph", dest="input_graph",
                  help="load graph (in JSON format)", metavar="FILE")
    parser.add_option("--pdb-id", dest="pdb_id",metavar="PDBID")
    parser.add_option("--filter-n-type", dest="filter_n_type",metavar="N_TYPE")
    parser.add_option("--filter-desc", dest="filter_desc",metavar="DESC")
    parser.add_option("--filter-full-desc", dest="filter_full_desc",metavar="DESC")

    (options, args)  = parser.parse_args()

    return (parser, options, args)

def main():
    (parser, options, args) = parse_args()
    assert options.input_graph is not None
    assert options.pdb_id is not None
    
    n_type_regexp = None
    desc_regexp = None
    full_desc_regexp = None
    if options.filter_n_type:
        n_type_regexp = re.compile('^'+options.filter_n_type+'$')
    if options.filter_desc:
        desc_regexp = re.compile('^'+options.filter_desc+'$')
    if options.filter_full_desc:
        full_desc_regexp = re.compile('^'+options.filter_full_desc+'$')
    
    g = GraphTool(options.input_graph)
    result = []
    for d_id in g.get_ids():
        full_id = options.pdb_id.upper() + ":"+d_id
        contacts = g.get_all_contacts_by_id(d_id,data=True)
        assert len(contacts)>0
        if options.filter_n_type:
            if not any([n_type_regexp.match(r.get('n_type','')) for r in contacts]):
                continue
        if options.filter_desc:
            if not any([desc_regexp.match(r.get('desc','')) for r in contacts]):
                continue
        if options.filter_full_desc:
            if not any([full_desc_regexp.match(r.get('full_desc','')) for r in contacts]):
                continue
        print full_id, contacts
        result.append(full_id)
    if options.output_json:
        save_json(options.output_json, result)

if __name__ == '__main__':
    main()

