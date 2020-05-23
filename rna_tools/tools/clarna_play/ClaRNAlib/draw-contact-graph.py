#!/usr/bin/env python
import sys
import os
import re
import math
import gzip
from optparse import OptionParser
from itertools import combinations
import shutil
import networkx as nx
from networkx.readwrite import json_graph

from utils import save_json

VARNA_HOME = "/Users/tomek/work/genesilico.programs/varna.gc"
if os.path.isdir("/home/twalen/prgs/varna.gc"):
    VARNA_HOME="/home/twalen/prgs/varna.gc"

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""draw the varna image for the given contact graph""")
    parser.add_option("-i", dest="input_graph", metavar="FILE")
    parser.add_option("-o", dest="output", metavar="FILE")
    parser.add_option("--resolution", dest="resolution", metavar="N", default=2)
    parser.add_option("--varna-home", dest="varna_home", metavar="DIR", default=VARNA_HOME)
    parser.add_option("--add-res-annotations", dest="add_res_annotations", action='store_true')
    parser.add_option("--add-idx-annotations", dest="add_idx_annotations", action='store_true')
    parser.add_option("--annotations-step", dest="annotations_step", metavar="N", default="10")
    parser.add_option("--ignore-edges", dest="ignore_edges", metavar="E1,E2,...")
    
    (options, args)  = parser.parse_args()
    options.annotations_step = int(options.annotations_step)
    return (parser, options, args)

def load_graph(fn):
    if fn is None:
        return None
    if not os.path.isfile(fn):
        print >>sys.stderr, "MISSING FILE %s" % fn
        return None
    if re.match(r".*\.gz$",fn):
        f = gzip.open(fn,"r")
    else:
        f = open(fn,"r")
    return json_graph.load(f)

def draw_contacts(graph, output_fn, options):
    nodes_dict = {}
    n = 0
    for node_id,d in graph.nodes(data=True):
        if len(d['resname'])>1:
            print "ignoring node %s, bad resname: %s" % (node_id, d['resname'])
            continue
        nodes_dict[node_id] = d
        n += 1
    nodes = sorted(nodes_dict.keys(), key=lambda x: "%c%09d" % (x[0], int(x[1:])))
    for i,_id in enumerate(nodes,start=1):
        nodes_dict[_id]['num']=i
    structure_dbn = "."*n
    sequence_dbn = "".join([nodes_dict[_id]['resname'] for _id in nodes])
    assert len(structure_dbn)==n
    assert len(sequence_dbn)==n
    aux_bps_list = []
    if options.ignore_edges:
        ignore_edges = options.ignore_edges.split(",")
    else:
        ignore_edges = []
    for v1,v2,d in sorted(graph.edges(data=True)):
        if d['type']=='contact' and d.get('reverse',False)==False:
            i1 = nodes_dict[v1]['num']
            i2 = nodes_dict[v2]['num']
            assert i1>=1 and i1<=n
            assert i2>=1 and i2<=n
            if "%s:%s"%(i1,i2) in ignore_edges or "%s:%s"%(v1,v2) in ignore_edges:
                continue
            if "?" in d['desc']:
                continue
            print "adding edge (%d,%d) - %s:%s desc=%s" % (i1,i2,v1,v2,d['desc'])
            if re.match(r'^[sSWH]{2}_(cis|tran)$', d['desc']):
                edge5 = d['desc'][0].upper()
                if edge5=='W':
                    edge5='WC'
                edge3 = d['desc'][1].upper()
                if edge3=='W':
                    edge3='WC'
                stericity = d['desc'][3:]
                if stericity=='tran':
                    stericity='trans'
                aux_bps_list.append("(%d,%d):edge5=%s,edge3=%s,stericity=%s,color=#4156C5" % (i1,i2,edge5,edge3,stericity))
            elif re.match(r'^[<>]{2}$', d['desc']):
                if d['desc']=='<<':
                    edge5 = ">>"
                    edge3 = ">>"
                    tmp_i = i1; i1=i2; i2=tmp_i
                else:
                    edge5 = d['desc']
                    edge3 = d['desc']
                aux_bps_list.append("(%d,%d):edge5=%s,edge3=%s,color=#BC8F8F" % (i1,i2,edge5,edge3))
            elif re.match(r'^[HWS]{1,2}_[0-9]', d['desc']):
                aux_bps_list.append("(%d,%d):edge5=B,edge3=P" % (i1,i2))

    aux_bps = ";".join(aux_bps_list)+";"
    
    resolution=int(options.resolution)
    
    jar = os.path.join(options.varna_home,"VARNA.jar")
    class_name = "fr.orsay.lri.varna.applications.VARNAcmd"
    cmd = "java -cp %(jar)s %(class_name)s -o %(output_fn)s -resolution %(resolution)d -warning true" % locals()
    cmd += " -structureDBN '%(structure_dbn)s'" % locals()
    cmd += " -sequenceDBN '%(sequence_dbn)s'" % locals()
    cmd += " -auxBPs '%(aux_bps)s'" % locals()
    if options.add_idx_annotations or options.add_res_annotations:
        annotations = []
        all_nodes = list(enumerate(nodes,start=1))
        for i,_id in [all_nodes[x] for x in range(0,len(nodes),options.annotations_step)]:
            l = ""
            if options.add_res_annotations:
                l += _id
            if options.add_idx_annotations:
                if l!="":
                    l += "-"
                l += str(i)
            annotations.append("%s:type=B,anchor=%d,size=8"%(l,i))
        cmd += " -annotations '%s'" % (";".join(annotations))

    print cmd
    os.system(cmd)

def main():
    (parser, options, args) = parse_args()
    graph = load_graph(options.input_graph)
    print "running varna"
    draw_contacts(graph, options.output, options)
    
if __name__ == '__main__':
    main()

