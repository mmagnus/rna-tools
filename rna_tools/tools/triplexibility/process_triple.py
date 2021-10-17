#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
0 = gray
1 = cHW AC does not exist	
2 = yellow
3 = green, exists

http://rna.bgsu.edu/triples/seq/CAC.html

for i in bucket/families/*; do python process_triple.py $i >> triple-db.csv ; done

"""
from __future__ import print_function
import argparse
import re

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args() # ['familes/triple_cWW_cHW.html'])
    v = args.verbose
    
    cols = """A C G U""".split()
    family = args.file.split('/')[-1].replace('.html', '').replace('triple_', '') # pages/AAA.html
    # collect all triples
    index = 0
    seqs = []
    """
    32 cWW_cHW_GAA to make a list that you can pick a triple seq
    like in seqs[0] to get the first one
    0 cWW_cHW_AAA
    1 cWW_cHW_AAC
    2 cWW_cHW_AAG
    3 cWW_cHW_AAU
    4 cWW_cHW_ACA
    5 cWW_cHW_ACC
    """
    for l in open(args.file):
        l = l.strip()
        if "<th class='narrow'" in l:
            for c in cols:
                l = l.replace("<th class='narrow'>", '').replace("</th>", '')
                if l:
                    triple = family + '_' + l + c #, end=' ')
                    if v: print(index, triple)
                    seqs.append(triple)
                    index += 1

    tds = open(args.file).read().split('<td')
    ntds = []
    data = []
    index = 0
    for td in tds[1:]:
        td = td.replace("class='wide small-font'>", '').replace('</td>', '')
        
        r = re.findall('<strong>(\d+)</strong> instance', td)
        instance = 0
        if r:
            if v: print(r[0])
            instance = int(r[0])

        clashes = 0
        r = re.findall('<strong>(\d+)</strong> clash', td)
        if r:
            if v: print(int(r[0]))
            clashes = int(r[0])

        near = 0
        r = re.findall("instances'>(\d+) near", td)
        if r:
            if v: print(int(r[0]))
            near = int(r[0])

        exists = 1
        if 'does not exist' in td:
            title = seqs[index]
            exists = 0

        r = re.findall("title='(.+)'>", td)
        title = ''
        if r:
            title = r[0].replace(' ', '_')
            if v: print(title)

        data.append([seqs[index], title, instance, clashes, near, exists])
        ntds.append(td)

        index += 1

    for f in data:
        print(','.join([str(x) for x in f]))
