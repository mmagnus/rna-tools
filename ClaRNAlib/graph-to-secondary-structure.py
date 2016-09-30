#!/usr/bin/env python
#
# Transforms contact graph into secondary structure (currently dot-bracket format)
# http://projects.binf.ku.dk/pgardner/bralibase/RNAformats.html
#

import sys
import os
import re
import gzip
import networkx as nx
from StringIO import StringIO
from optparse import OptionParser
# biopython
from Bio import PDB
#
from utils import write_file, load_graph, load_pdb

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="transform contact graph into secondary stucture")
    parser.add_option("-o", "--output", dest="output",
                  help="save output to FILE", metavar="FILE")
    parser.add_option("-i", "--input", dest="input",
                  help="read input graph from FILE", metavar="FILE")
    parser.add_option("--name", dest="name",
                  help="name of the structure", metavar="NAME")
    parser.add_option("--pdb", dest="pdb",
                  help="read PDB from FILE", metavar="FILE")
    (options, args)  = parser.parse_args()
    return (parser, options, args)

def get_nodes_ordering(graph, options):
    if options.pdb:
        res = []
        pdb = load_pdb(options.pdb)
        graph_nodes = set(graph.nodes())
        for c in pdb.get_chains():
            for r in c:
                node_id = c.id+str(r.get_id()[1])+r.get_id()[2].strip()
                if node_id in graph_nodes:
                    res.append(node_id)
        return res
    else:
        def nodes_ordering(x):
            m = re.match("^(.)([0-9]+)(.*)$", x)
            if m:
                return (m.group(1), int(m.group(2)), m.group(3))
            else:
                return (x[0],x[1:],"")
        return sorted(graph.nodes(), key=nodes_ordering)

def get_sequence(graph, nodes):
    res = []
    return "".join(graph.node[n].get('resname','?') for n in nodes)
    return ""

def get_structure_name(options):
    if options.pdb:
        return re.sub("(.pdb)?(.gz)?$", "", os.path.basename(options.pdb))
    else:
        return re.sub("(.json)?(.gz)?$", "", os.path.basename(options.input))

def gen_secondary_structure(g, options):

    name = get_structure_name(options)
    nodes = get_nodes_ordering(g, options)
    nodes2idx = dict([(n, i) for i, n in enumerate(nodes)])

    sequence = get_sequence(g, nodes)

    gg = nx.DiGraph()
    canonical_pairs = ('CG','GC','AU','UA')
    for u,v,data in g.edges(data=True):
        # print u,v,data
        if data.get('desc') == 'WW_cis' and data.get('n_type') in canonical_pairs:
            gg.add_edge(u, v)
    brackets = ['.' for n in nodes]
    for v1, v2 in nx.algorithms.maximal_matching(gg):
        i1 = nodes2idx[v1]
        i2 = nodes2idx[v2]
        brackets[min(i1, i2)] = '('
        brackets[max(i1, i2)] = ')'

    res = ""
    res += "> "+name+"\n"
    res += sequence+"\n"
    res += "".join(brackets)+"\n"
    return res

def main():
    (parser, options, args) = parse_args()
    if options.input and options.output:
        g = load_graph(options.input)
        res = gen_secondary_structure(g, options)
        write_file(options.output, res)
    else:
        print "specify input and output"
        parser.print_help()
        exit(1)

if __name__ == '__main__':
    main()

