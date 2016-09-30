#!/usr/bin/python
import sys
import os
import re
import math
import random
from optparse import OptionParser
import simplejson as json

from numpy import array
from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from StringIO import StringIO
import scipy as scipy
import numpy as np
import scipy.cluster.hierarchy as sch

from utils import load_json, save_file

def cluster_motifs(options):
    rmsd_tolerance = float(options.rmsd_tolerance)

    motifs = pdb_files(options.input_dir)
    neg_motifs = pdb_files(options.input_dir,include_positive=False,include_negative=True)
    parser = PDB.PDBParser()

    clusters = {}

    for m in motifs:
        mm = os.path.basename(m)
        key1 = mm[0:2]
        key2 = re.sub(r"^(..)_desc_(.*)_from.*$",r"\2",mm)
        # print m, key1, key2
        s = parser.get_structure("c", m)
        s = simplify_structure_old(s, key1)
        if s is None:
            print "ignoring: %s" % m
            continue
        if not clusters.has_key(key1):
            clusters[key1] = []
        similar = None
        for (i,(c,count,desc,min_neg_distance)) in enumerate(clusters[key1]):
            if options.atoms:
                d = new_res_distance(c[0],c[1], s[0], s[1])
            else:
                d = res_distance(c[0],c[1], s[0], s[1])
            if d < rmsd_tolerance:
                similar = i
                break
        if similar is not None:
            clusters[key1][similar][1] += 1
            if not clusters[key1][similar][2].has_key(key2):
                clusters[key1][similar][2][key2] = 0
            clusters[key1][similar][2][key2] += 1
            print "new element in cluster %d, new size: %d" % (similar, clusters[key1][similar][1])
        else:
            print "new cluster %s" % m
            clusters[key1].append([s,1,{key2:1},1000])

    neg_distances = {}
    for m in neg_motifs:
        mm = os.path.basename(m)
        key1 = mm[0:2]
        if not clusters.has_key(key1):
            continue
        if not neg_distances.has_key(key1):
            neg_distances[key1] = [1000]*len(clusters[key1])
        s = parser.get_structure("c", m)
        s = simplify_structure_old(s, key1)
        min_i = None
        min_d = None
        for (i,(c,count,desc,min_neg_distance)) in enumerate(clusters[key1]):
            if options.atoms:
                d = new_res_distance(c[0],c[1], s[0], s[1])
            else:
                d = res_distance(c[0],c[1], s[0], s[1])
            if d < neg_distances[key1][i]:
                print "min. negative distance for cluster %s/%d has been decrased to %.5f" % (key1,i, d)
                neg_distances[key1][i] = d
                clusters[key1][i][3] = d

    total = 0
    for k in clusters.keys():
        total += len(clusters[k])
        print "%s : %d clusters" % (k, len(clusters[k]))
        for (i,(c,size,desc,neg_distance)) in enumerate(clusters[k]):
            print " #%04d, size=%d, min negative distance=%.5f" % (i, size, neg_distances[k][i])
            for d,size2 in desc.items():
                print "   %s: %d" % (d,size2)
    print "total clusters: %d" % total
    if options.output:
        f = open(options.output,"w")
        json.dump(clusters, f)
        f.close()

def atoms_list(atoms_dict, selected_atoms):
    return [atoms_dict.get(x) for x in selected_atoms]

def make_doublet(r1,r2):
    SELECTED_ATOMS = ['C2','C4','C6',"C1'"]
    vec = atoms_list(r1['atoms'],SELECTED_ATOMS) + atoms_list(r2['atoms'],SELECTED_ATOMS)
    if None in vec:
      return None
    vec2 = atoms_list(r2['atoms'],SELECTED_ATOMS) + atoms_list(r1['atoms'],SELECTED_ATOMS)
    return {
        'id': r1['id']+':'+r2['id'],
        'n_type': r1['resname']+r2['resname'],
        'vec': array(vec,'f'),
        'vec2': array(vec2,'f'),
    }

def load_doublets(fn,options):
    doublets = load_json(fn)
    if options.limit:
        print "processing only first %s doublets" % options.limit
        doublets = doublets[0:int(options.limit)]
    doublets_dict = {}
    by_pdb = {}
    for x in doublets:
        xx = x.split(":")
        if not by_pdb.has_key(xx[0]):
            by_pdb[xx[0]] = set()
        by_pdb[xx[0]].add(x)
    print "loaded %s doublets from %d pdb files" % (len(doublets), len(by_pdb.keys()))
    for pdb in sorted(by_pdb.keys()):
        r_fn = "%s/%s_rna.residues.json.gz" % (options.residues_dir,pdb.lower())
        assert os.path.isfile(r_fn)
        residues = load_json(r_fn)
        print " - %s, residues: %s, doublets: %s" % (pdb, len(residues), len(by_pdb[pdb]))
        for doublet in by_pdb[pdb]:
            xx = doublet.split(":")
            if not residues.has_key(xx[1]) or not residues.has_key(xx[2]):
                print >> sys.stderr, "MISSING doublet %s" % doublet
                continue
            assert residues.has_key(xx[1])
            assert residues.has_key(xx[2])
            d = make_doublet(residues[xx[1]], residues[xx[2]])
            if d is None:
                print >> sys.stderr, "WRONG doublet %s" % doublet
                continue
            doublets_dict[doublet] = d
    return doublets, doublets_dict

def doublets_dist(d1, d2):
    sup = SVDSuperimposer()
    sup.set(d1['vec'], d2['vec'])
    sup.run()
    rms1 = sup.get_rms()
    sup.set(d1['vec'], d2['vec2'])
    sup.run()
    rms2 = sup.get_rms()
    return min(rms1, rms2)

def group_doublets(doublets, doublets_dict, max_rmsd):
    new_doublets = []
    new_doublets_elems = []
    new_weights = []
    for i,d in enumerate(doublets):
        if not doublets_dict.has_key(d):
            continue
        if i%100==0:
            print "i=%d, new_doublets=%d" % (i,len(new_doublets))
            sys.stdout.flush()
        k = None
        for j,dd in enumerate(new_doublets):
            if doublets_dist(doublets_dict[d], doublets_dict[dd]) <= max_rmsd:
                k = j
                break
        if k is None:
            new_doublets.append(d)
            new_doublets_elems.append([d])
            new_weights.append(1)
        else:
            new_weights[k] += 1
            new_doublets_elems[k].append(d)
    return new_doublets, new_doublets_elems, new_weights

def compute_dist_matrix(doublets, doublets_dict):
    n = len(doublets)
    res = [ [0]*n for i in xrange(n)]
    for i in xrange(n-1):
        print "computing dist matrix i=%s n=%s" % (i,n)
        sys.stdout.flush()
        for j in xrange(i+1,n):
            d = doublets_dist(doublets_dict[doublets[i]], doublets_dict[doublets[j]])
            res[i][j] = d
            res[j][i] = d
    return res

def cluster_doublets(options):
    doublets, doublets_dict = load_doublets(options.input_json, options)
    (new_doublets, new_doublets_elems, new_weights) = group_doublets(doublets, doublets_dict, options.rmsd_tolerance)
    print "after pruning %s" % len(new_doublets)

    print "computing dist matrix"
    sys.stdout.flush()
    m = compute_dist_matrix(new_doublets, doublets_dict)
    data = np.array(m, 'f')
    print "done"
    sys.stdout.flush()

    condensed = scipy.spatial.distance.squareform(data)

    fig = plt.figure(figsize=(8,8))

    print "computing clustering"
    Y = sch.linkage(condensed, method='complete', metric='euclidean')
    sys.stdout.flush()
    clusters_info = sch.fcluster(Y, t=options.max_cluster_rmsd, criterion='distance')

    ax1 = fig.add_axes([0.07,0.1,0.1,0.7])
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])

    ax2 = fig.add_axes([0.2,0.83,0.7,0.17])
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])

    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    data = data[idx1,:]
    data = data[:,idx2]

    axmatrix = fig.add_axes([0.2,0.1,0.7,0.7])
    im = axmatrix.matshow(data, aspect='auto', origin='lower', cmap=plt.cm.Blues)

    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    axcolor = fig.add_axes([0.93,0.1,0.02,0.7])
    plt.colorbar(im, cax=axcolor)

    f = StringIO()
    fig.savefig(f, format="png")
    plt.close(fig)
    res = f.getvalue()

    if options.save_png:
        print "saving: %s" % options.save_png
        save_file(options.save_png, res)
    if options.save_json:
        print "saving: %s" % options.save_json
        json.dump(
            {
                'doublets': new_doublets,
                'doublets_elems': new_doublets_elems,
                'idx': idx1,
                'clusters_info': clusters_info.tolist()
            },
            open(options.save_json,"w"),
            indent=2
        )


def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""cluster doublets, used for clustering unclassified doublets""")
    parser.add_option("-i", "--input-json", dest="input_json",
                  help="read input pairs from JSON file", metavar="JSON")
    parser.add_option("--limit", dest="limit",
                  help="process only first X doublets", metavar="X")
    parser.add_option("--residues-dir", dest="residues_dir",
                  help="the directory with the residues", metavar="DIR")
    parser.add_option("--rmsd-tolerance", dest="rmsd_tolerance",
                  help="rmsd tolerance for prunning input set", metavar="RMSD",
                  default="0.5")
    parser.add_option("--max-cluster-rmsd", dest="max_cluster_rmsd",
                  help="maximal rmsd in single cluster", metavar="RMSD",
                  default="1.0")
    parser.add_option("--save-png", dest="save_png",
                  help="save clustering results to PNG file", metavar="FILE")
    parser.add_option("--save-json", dest="save_json",
                  help="save clustering results to JSON file", metavar="JSON")

    (options, args)  = parser.parse_args()
    options.rmsd_tolerance = float(options.rmsd_tolerance)
    options.max_cluster_rmsd = float(options.max_cluster_rmsd)
    return (parser, options, args)

def main():
    (parser, options, args) = parse_args()
    cluster_doublets(options)

if __name__ == '__main__':
    main()


