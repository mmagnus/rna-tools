#!/usr/bin/env python
#
import re
import sys
import itertools
import multiprocessing
import random
from optparse import OptionParser

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scipy as scipy
import numpy as np
import scipy.cluster.hierarchy as sch
from utils import *
from distances import rmsd_distance, doublet_params_dict, normalize_points, Doublet, Residue
from progressbar import ProgressBar, Percentage, Bar, ETA

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""some experiments with classifier generation""")
    parser.add_option("--doublet-id", dest="doublet_id",
                  help="show params for doublet", metavar="PDBID:R1:R2")
    parser.add_option("--input-pdb", dest="input_pdb",
                  help="read doublets from PDB file", metavar="FILE")
    parser.add_option("--input-json", dest="input_json",
                  help="read doublets from json file", metavar="FILE")
    parser.add_option("--only-keys", dest="only_keys",
                  help="process only following keys from input file", metavar="KEYS")
    parser.add_option("--input-group-info", dest="input_group_info",
                  help="read doublets from group info file", metavar="FILE")
    parser.add_option("--show-histograms", dest="show_histograms", action='store_true',
                  help="show histograms", default=False)
    parser.add_option("--show-2d-histograms", dest="show_2d_histograms", action='store_true',
                  help="show 2D histograms", default=False)
    parser.add_option("--data-dir", dest="data_dir",
                  help="directory with data", metavar="DIR", default="gc-data")
    parser.add_option("--params-type", dest="params_type",
                  help="parameters type", metavar="T", default="bp")
    parser.add_option("--eval-for", dest="eval_for",
                  help="process only keys for given DESC", metavar="DESC")
    parser.add_option("--gen-pdb-for", dest="gen_pdb_for",
                  help="gen pdb file")
    parser.add_option("-o","--output", dest="output",
                  help="save output to file", metavar="FILE")
    parser.add_option("--limit", dest="limit",
                  help="limit number of elements", metavar="N")
    (options, args)  = parser.parse_args()
    if options.gen_pdb_for:
        if not options.input_group_info:
            prefix = "unk"
            if "BPh" in options.gen_pdb_for:
                prefix = "base-phosphate"
            elif "BR" in options.gen_pdb_for:
                prefix = "base-ribose"
            elif "<" in options.gen_pdb_for or ">" in options.gen_pdb_for:
                prefix = "stacking"
            else:
                prefix = "bp"
            options.input_group_info = os.path.join(options.data_dir,"cl",prefix+"_"+options.gen_pdb_for.replace("/","_")+".json.gz")
            options.n_type = options.gen_pdb_for.split("/")[1]
    elif options.eval_for:
        if not options.only_keys:
            options.only_keys = ",".join(["evaluation/"+k+"/"+options.eval_for for k in ('ref-ok','ref-ok-fuzzy','ref-undetected','ref-diff','prev-undetected','ref-all')])
        if not options.input_json:
            options.input_json = FileNamesObject.eval_fn(setname="training",t="eval1")
    if options.limit is not None:
        options.limit = int(options.limit)
    return (parser, options, args)

def compute_params(doublet_lists, labels, options):
    assert len(doublet_lists) == len(labels)
    if options.limit:
        for label,x in zip(labels,doublet_lists):
            if len(x)>options.limit:
                print "reducing number of elements of set %s from %d to %d" % (label,len(x),options.limit)
        doublet_lists = [x[0:options.limit] for x in doublet_lists]
    l_count = len(doublet_lists)
    all_ids = set(sum(doublet_lists,[]))
    r_atoms = ['P',"C1'","O2'","O3'","O4'"]
    r_atoms += ["OP1","OP2","O5'","NEXT:O3'"]
    # base atoms
    r_atoms += ['N9', 'C4', 'N3', 'N1', 'C6', 'N6', 'C8', 'C5', 'C2', 'N7']
    r_atoms += ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C6', 'C5']
    r_atoms += ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C6', 'C5']
    r_atoms += ['N9', 'C4', 'N3', 'N1', 'C6', 'O6', 'C8', 'C5', 'C2', 'N7', 'N2']
    r_atoms = set(r_atoms)
    
    if options.params_type in ['base-ribose']:
        r_atoms += ["C4'","C3'","C2'","C1'"]
    dd = DoubletsDict(options.data_dir, reduced_atoms=r_atoms)
    dd.load_pdb_files(all_ids, verbose=True)
    
    res = {}
    stats = [{} for i in xrange(l_count)]
    total = [0 for i in xrange(l_count)]
    for i,d_list in enumerate(doublet_lists):
        for d_id in d_list:
            print d_id
            d = dd.get(d_id)
            if d is None or d[0] is None or d[1] is None:
                print "INVALID doublet! %s" % d_id
                continue
            n_type = dd.get_n_type(d_id)
            # p = doublet_params_dict(d, n_type, options.params_type)
            p = Doublet(d_id,Residue("A1",n_type[0],d[0]),Residue("B1",n_type[1],d[1])).get_all_params()
            # p['dist_sum']=p['min_dist']*4+p['dist']
            print i, labels[i], d_id, n_type, [(k,p[k]) for k in sorted(p.keys())]
            for k,v in p.items():
                if not res.has_key(k):
                    res[k] = [[] for j in xrange(l_count)]
                res[k][i].append(v)
            total[i] += 1
            key = " ".join([k+"="+str(v) for k,v in sorted(p.items()) if re.match('^i_',k) or k=='oxygens_count' or re.match('^strand_orient',k)])
            if not stats[i].has_key(key):
                stats[i][key] = []
            stats[i][key].append(d_id)
    for i in xrange(l_count):
        for k in sorted(stats[i].keys(), key=lambda x: len(stats[i][x]), reverse=True):
            print "STATS[%s] %s count=%d pr=%.2f%% samples=%s" % (labels[i], k,len(stats[i][k]), float(100*len(stats[i][k]))/total[i], stats[i][k][0:10])
    if options.show_2d_histograms:
        for k1,k2 in [('nn_ang_norm','n12cc_ang'),('rot_ang','min_dist'),('dist','min_dist'),('n2_z','nn_ang_norm')]:
            l = len(res[k1])
            for i in xrange(l):
                data_x = []
                data_y = []
                for j in xrange(len(res[k1][i])):
                    x = res[k1][i][j]
                    y = res[k2][i][j]
                    if x is not None and y is not None:
                        data_x.append(x)
                        data_y.append(y)
                        print labels[i], x, y
                H, xedges, yedges = np.histogram2d(data_x, data_y, bins=(20,20), range=[[min(data_x),max(data_x)], [min(data_y),max(data_y)]], normed=False)
                # czy tu sa dobre wartosci?
                plt.pcolor(xedges, yedges, H, cmap=matplotlib.cm.get_cmap('OrRd'))
                # extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                # plt.imshow(H.transpose(), extent=extent, interpolation='nearest', origin="lower", cmap=matplotlib.cm.get_cmap('OrRd'))
                plt.colorbar()
                plt.xlabel(k1)
                plt.ylabel(k2)
                plt.title("(%s,%s) - %s" % (k1,k2,labels[i]))
                plt.show()
            # draw all series on single plot
            if l>1:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                colors = ('r', 'g', 'b', 'k')
                for i in xrange(l):
                    data_x = []
                    data_y = []
                    for j in xrange(len(res[k1][i])):
                        x = res[k1][i][j]
                        y = res[k2][i][j]
                        if x is not None and y is not None:
                            data_x.append(x)
                            data_y.append(y)
                    ax.scatter(data_x, data_y, marker='o', c=colors[i%len(colors)], label=labels[i])
                plt.legend()
                ax.set_xlabel(k1)
                ax.set_ylabel(k2)
                plt.title("(%s,%s)" % (k1,k2))
                plt.show()
        

    if options.show_histograms:
        for k in sorted(res.keys()):
            data = [[x for x in row if x is not None] for row in res[k]]
            sum_data = sum(data,[])
            if len(sum_data)==0:
                print "empty data set"
                continue
            if k in ['oxygens','br_info','ph_info']:
                continue
            
            min_v = min(sum_data)
            max_v = max(sum_data)
            print "k", k, "min_v", min_v, "max_v", max_v, "data"
        
            fig = plt.figure()
            n, bins, patches = plt.hist(data, 25, range=(min_v,max_v), normed=0, histtype='bar',
                            label=labels)
            plt.legend()
            plt.title(k)
            plt.show()

def points2model(model_num,points,n_type):
    num = 0
    model = PDB.Model.Model(model_num)
    for i,p in enumerate(points):
        chain = PDB.Chain.Chain(chr(ord('A')+i))
        num += 1
        residue = PDB.Residue.Residue((' ',i,' '),n_type[i],num)
        for j,(k,v) in enumerate(p.items()):
            num += 1
            kk = k.replace("NEXT:","")
            element = k.replace("NEXT:","").replace("1","").replace("2","")[0]
            # print "serial: %s, k=%s, element=%s" % (num,k,element)
            atom = PDB.Atom.Atom(kk,v,1,1," ",fullname=kk,serial_number=num,element=element)
            residue.add(atom)
        chain.add(residue)
        model.add(chain)
    return model

def gen_pdb(doublets,doublets_unclassified,options):
    from distances import NORMALIZED_BASE
    from Bio.SVDSuperimposer import SVDSuperimposer

    n_type0 = options.n_type[0]

    SUP_ATOMS = {
      'A': ['N9', 'C4', 'N3', 'N1', 'C6', 'N6', 'C8', 'C5', 'C2', 'N7'],
      'C': ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C6', 'C5'],
      'U': ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C6', 'C5'],
      'G': ['N9', 'C4', 'N3', 'N1', 'C6', 'O6', 'C8', 'C5', 'C2', 'N7', 'N2'],
    }


    r_atoms = ['N1','N2','N3','N4','N6','N7','P']
    r_atoms += ['C2','C4','C5','C6','C8']
    r_atoms += ['O2']
    r_atoms += ["O2'","O3'","O4'"]
    r_atoms += ["OP1","OP2","O5'","NEXT:O3'"]
    r_atoms += ['N1','C6','O6','C5','C4','N3','C2','N2','N7','C8','N9']

    #  'A': ['N9', 'C4', 'N3', 'N1', 'C6', 'N6', 'C8', 'C5', 'C2', 'N7'],
    #  'C': ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C6', 'C5'],
    #  'U': ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C6', 'C5'],
    #  'G': ['N9', 'C4', 'N3', 'N1', 'C6', 'O6', 'C8', 'C5', 'C2', 'N7', 'N2'],

    dd = DoubletsDict(options.data_dir, reduced_atoms=r_atoms)
    dd.load_pdb_files(doublets+doublets_unclassified, verbose=True)
    
    for d_set,suffix in (doublets,""),(doublets_unclassified,"-uncl"):
        s = PDB.Structure.Structure("superimpose")
        num = 0
        model = points2model(1,(NORMALIZED_BASE[n_type0],{}),options.n_type)
        s.add(model)
        for d_num,d_id in enumerate(d_set,start=2):
            print d_num,d_id
            new_p = normalize_points(dd.get(d_id), n_type0)
            for k,v in new_p[0].items():
                if k not in SUP_ATOMS[n_type0]:
                    del new_p[0][k]
            for k,v in new_p[1].items():
                if k not in ['OP1','OP2','P',"O5'","NEXT:O3'"]:
                    del new_p[1][k]
            model = points2model(d_num,new_p,options.n_type)
            s.add(model)
        out_fn = options.output.replace(".pdb","")+suffix+".pdb"
        save_pdb(out_fn,s)

########################

def main():
    (parser,options,_args) = parse_args()
    
    doublet_lists = []
    labels = []
    if options.doublet_id:
        doublet_lists.append([options.doublet_id])
        labels.append("d:%s" % options.doublet_id)
    elif options.gen_pdb_for:
        json = load_json(options.input_group_info)
        assert isinstance(json,dict)
        assert len(json.keys())==1
        v = []
        vu = []
        for group_info in json.values():
            v = group_info['all_doublets']
            vu = list(set([x for row in group_info['neigh_unclassified'] for x,d in row]))
        gen_pdb(v,vu,options)
    elif options.input_group_info:
        json = load_json(options.input_group_info)
        assert isinstance(json,dict)
        assert len(json.keys())==1
        for group_info in json.values():
            assert isinstance(group_info,dict)
            for k in ('all_doublets', 'neigh_unclassified', 'neigh_other'):
                v = group_info[k]
                if k=='all_doublets':
                    doublet_lists.append(v)
                    labels.append('reference')
                else:
                    vv = list(set([did for row in v for (did,dist) in row]))
                    doublet_lists.append(vv)
                    labels.append(k.split("_")[1])
    elif options.input_json:
        only_keys = None
        if options.only_keys:
            only_keys = options.only_keys.split(",")

        for fn in options.input_json.split(","):
            if fn=='':
                continue
            json = load_json(fn)
            if isinstance(json,dict):
                print "DICT!"
                keys = json.keys()
                if only_keys is not None:
                    if len(only_keys)==1:
                        regexp = re.compile('^'+only_keys[0]+'$')
                        keys = [k for k in json.keys() if regexp.match(k)]
                    else:
                        keys = only_keys
                for k in keys:
                    if not json.has_key(k):
                        continue
                    v = json[k]
                    assert all([isinstance(did,str) for did in v])==True
                    print k
                    doublet_lists.append(v)
                    labels.append(k)
                if only_keys is not None and len(only_keys)==1:
                    doublet_lists = [sum(doublet_lists,[])]
                    labels = [only_keys[0]]
            elif isinstance(json,list):
                print "LIST!"
                assert all([isinstance(did,str) for did in json])==True
                doublet_lists.append(json)
                labels.append(os.path.basename(fn))
            else:
                raise Exception("Unknown format of JSON file")
    elif options.input_pdb:
        structure = load_pdb(options.input_pdb)
        residues = [r for r in structure.get_residues()]
        assert len(residues)==2
        n_type = residues[0].resname.strip()+residues[1].resname.strip()
        print doublet_params_dict((simplify_residue(residues[0]),simplify_residue(residues[1])),n_type,options.params_type)
    
    compute_params(doublet_lists, labels, options)
    
if __name__=="__main__":
    main()