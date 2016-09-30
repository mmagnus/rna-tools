#!/usr/bin/env python
"""make dictionary of reverse descriptions"""
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

from utils import *

CMD_PYTHON="/usr/bin/python"
if os.path.isfile("/opt/local/bin/python2.7"):
    CMD_PYTHON = "/opt/local/bin/python2.7"

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="create JSON dictionary for reverse descriptions")
    parser.add_option("-o", "--output-json", dest="output_json",
                  help="save dictionary to JSON file", metavar="DIR")
    parser.add_option("-i", "--input", dest="input",
                  help="process given PDB files", metavar="FILE")
    (options, args)  = parser.parse_args()

    return (parser, options, args)


def handle(options):
    inp_fn = options.input
    res_map = inp_fn.replace(".pdb.gz",".json.gz")
    res = {}
    (_fd,tmp_fn) = tempfile.mkstemp(".make-rev-dict.json")
    for prg,prg_name in (('RV','rnaview'),('MC','mc-annotate'),('MO','moderna'),('FR','fr3d')):
        for t in ('bp','stackings'):
            if t=='bp' and prg=='MO':
                continue
            cmd = CMD_PYTHON + " ./extract-contacts.py --use-%(prg_name)s -i '%(inp_fn)s' --no-pdb --output-json=%(tmp_fn)s" % locals()
            if t=='stackings':
                cmd += " --extract-stackings"
            if prg=='FR':
                cmd += " --dont-normalize --residue-mapping=%(res_map)s" % locals()

            os.system(cmd + "> /dev/null")
            if os.path.exists(tmp_fn):
                doublets = load_json(tmp_fn)
            else:
                print "missing doublets dict!"
                doublets = []
            os.unlink(tmp_fn)

            if prg!='FR':
                os.system(cmd+" --reverse-pdb > /dev/null")
                rev_doublets = load_json(tmp_fn)
                os.unlink(tmp_fn)
            else:
                rev_doublets = []

            descriptions = {}
            for d in doublets+rev_doublets:
                id1 = d['resseq_1']
                id2 = d['resseq_2']
                if id1==id2:
                    continue
                desc = d['desc']
                key = (id1,id2)
                if not descriptions.has_key(key):
                    descriptions[key] = desc
            for (id1,id2),desc in descriptions.items():
                rev_key = (id2,id1)
                if descriptions.has_key(rev_key):
                    rev_desc = descriptions[rev_key]

                    key1 = "reverse/%s/%s/%s" % (prg, desc, rev_desc)
                    key2 = "reverse/%s/%s/%s" % (prg, rev_desc, desc)
                    if not res.has_key(key1):
                        res[key1] = 0
                    if not res.has_key(key2):
                        res[key2] = 0
                    res[key1] += 1
                    res[key2] += 1

    save_json(options.output_json, res)

def main():
    (parser, options, args) = parse_args()
    if not options.input:
        print "select input file"
        parser.print_help()
        exit(1)
    if not options.output_json:
        print "select output"
        parser.print_help()
        exit(1)
    handle(options)

if __name__ == '__main__':
    main()

