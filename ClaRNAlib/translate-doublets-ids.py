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
    parser = OptionParser(description="""create JSON dictionary with doublets
""")
    parser.add_option("-o", "--output-json", dest="output_json",
                  help="save result to JSON file", metavar="JSON")
    parser.add_option("-i", "--input-json", dest="input_json",
                  help="read input from JSON", metavar="JSON")
    parser.add_option("--input-str", dest="input_str",
                  help="translate doublets from string", metavar="STR")
    parser.add_option("--data-dir", dest="data_dir",
                  help="directory with data", metavar="DIR", default="gc-data")                  
    (options, args)  = parser.parse_args()

    return (parser, options, args)

def main():
    (parser, options, args) = parse_args()
    
    if options.input_json:
        input = load_json(options.input_json)
    elif options.input_str:
        input = [x.strip() for x in options.input_str.split(",")]
    else:
        print "select input"
        parser.print_help()
        exit(1)
        
    dd = DoubletsDict(options.data_dir,reduced_atoms=[])

    if isinstance(input,list):
        dd.load_pdb_files(input,verbose=True)
        res = [dd.get_org_id(did) for did in input]
    elif isinstance(input,dict):
        dd.load_pdb_files(set(sum(input.values(),[])),verbose=True)
        res = {}
        for k,v in input.items():
            res[k] = [dd.get_org_id(did) for did in v]
    else:
        print "unsupported input format"
        exit(1)
    
    if options.output_json:
        save_json(options.output_json, res)
    else:
        print res

if __name__ == '__main__':
    main()

