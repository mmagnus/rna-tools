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
    parser = OptionParser(description="""helper tool for minipulating JSON files""")
    parser.add_option("-i", "--input-json", dest="input_json",
                  help="input json file", metavar="FILE")
    parser.add_option("-o", "--output-json", dest="output_json",
                  help="output json file", metavar="FILE")
    parser.add_option("--output-indent", dest="output_indent",
                  help="output json indent", metavar="N")
    parser.add_option("--apply-map-op", dest="map_op",
                  help="apply map operation: (copy|copy-values|len|sort|print-values|print-keys|print-lengths)", metavar="OP")
    parser.add_option("--only-keys", dest="only_keys",
                  help="process only following keys from input file", metavar="KEYS")
    parser.add_option("--rename-keys", dest="rename_keys",
                  help="rename keys", metavar="IS,SHOULD_BE")
    parser.add_option("--skip-empty", dest="skip_empty", action='store_true',
                  default=False, help="skip creation of empty files")
    (options, args)  = parser.parse_args()

    return (parser, options, args)

def main():
    (parser, options, args) = parse_args()
    
    assert options.input_json is not None
    assert options.map_op in ['len','print-values','print-keys','print-lengths','copy','copy-values']
    
    json = load_json(options.input_json)
    if options.only_keys:
        only_keys = set(options.only_keys.split(","))
        for key in json.keys():
            if not any([re.match(o,key) for o in only_keys]):
                del json[key]
    if options.rename_keys:
        s_is,s_should_be = options.rename_keys.split(",")
        json = dict([(re.sub(s_is,s_should_be,k),v) for k,v in json.items()])

    output_indent = None
    if options.output_indent:
        output_indent = int(options.output_indent)

    if options.skip_empty:
        if json is None or len(json)==0:
            print "skipping empty file"
            return
    
    if options.map_op=='print-keys':
        for key in sorted(json.keys()):
            print key
    elif options.map_op=='print-values':
        for key in sorted(json.keys()):
            print "%s=%s" % (key,json[key])
    elif options.map_op=='print-lengths':
        for key in sorted(json.keys()):
            print "%s=%s" % (key,len(json[key]))
    elif options.map_op=='len':
        assert options.output_json is not None
        json = dict([(k,len(v)) for k,v in json.items()])
        save_json(options.output_json, json, indent=output_indent)
    elif options.map_op=='copy':
        assert options.output_json is not None
        save_json(options.output_json, json, indent=output_indent)
    elif options.map_op=='copy-values':
        assert options.output_json is not None
        save_json(options.output_json, sum(json.values(),[]), indent=output_indent)
    else:
        raise Exception("Unknown map_op")

if __name__ == '__main__':
    main()

