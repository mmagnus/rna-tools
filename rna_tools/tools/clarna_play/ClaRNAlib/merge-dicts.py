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
from itertools import combinations,product

from utils import *

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="merge multiple JSON dictionaries to single JSON dictionary")
    parser.add_option("-o", "--output", dest="output_json",
                  help="save dictionary to JSON file", metavar="FILE")
    parser.add_option("--input", dest="input",
                  help="comma separated list of files", metavar="FILE1,FILE2,...")
    parser.add_option("--input-from-file", dest="input_from_file",
                  help="read input file names from file", metavar="FILE")
    parser.add_option("--only-keys", dest="only_keys",
                  help="process only following keys from input file", metavar="KEYS")
    (options, args)  = parser.parse_args()

    return (parser, options, args)

def handle_merge(options,args):
    res = None

    files = []
    if options.input:
        files += options.input.split(",")
    if options.input_from_file:
        files += read_file(options.input_from_file).split("\n")
    if len(args)>0:
        files += args
        
    only_keys = None
    if options.only_keys:
        only_keys = set(options.only_keys.split(","))

    for fn in files:
        if fn=="":
            continue
        print "loading %s" % fn
        sys.stdout.flush()
        r = load_json(fn)
        if isinstance(r,list):
            if res is None:
                res = []
            res += r
        elif isinstance(r,dict):
            if res is None:
                res = {}
            for k,v in r.items():
                if only_keys is not None:
                    if not any([re.match('^'+o+'$',k) for o in only_keys]):
                        continue
                if not res.has_key(k):
                    res[k] = v
                else:
                    res[k] += v
        else:
            raise Exception("Unsupported content type")
    if res is None:
        res = {}
    save_json(options.output_json, res)

def main():
    (parser, options, args) = parse_args()
    if not options.output_json:
        print "select output"
        parser.print_help()
        exit(1)
    handle_merge(options, args)

if __name__ == '__main__':
    main()

