#!/usr/bin/env python
#
# script for computing the classifier
#
import re
import sys
import itertools
import multiprocessing
import random
from optparse import OptionParser
import scipy as scipy
import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.cluster.vq import kmeans, whiten
from scipy.spatial import KDTree
from utils import *
from distances import rmsd_distance, bp_distance, bph_distance, doublet_params_dict, \
    normalize_points, vector_length, center_vector, \
    expected_strand_orient, expected_stack_orient, \
    PH_OXYGENS, BR_OXYGENS
from progressbar import ProgressBar, Percentage, Bar, ETA

OUTPUT_DIR = "gc-data/comp_cl"
N_TYPES = [n1+n2 for n1,n2 in itertools.product(["A","C","U","G"],repeat=2)]

R_ATOMS = ['P',"C1'","O2'","O3'","O4'"]
R_ATOMS += ["OP1","OP2","O5'","NEXT:O3'"]
R_ATOMS += ["C4'","C3'","C2'","C1'","O4'"] # for ribose-center vector
# base atoms
R_ATOMS += ['N9', 'C4', 'N3', 'N1', 'C6', 'N6', 'C8', 'C5', 'C2', 'N7']
R_ATOMS += ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C6', 'C5']
R_ATOMS += ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C6', 'C5']
R_ATOMS += ['N9', 'C4', 'N3', 'N1', 'C6', 'O6', 'C8', 'C5', 'C2', 'N7', 'N2']

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""some experiments with classifier generation""")
    parser.add_option("--groups", dest="groups",
                  help="groups filename", metavar="FILE", default=FileNamesObject.groups_fn(setname="training",reduced=True))
    parser.add_option("--data-dir", dest="data_dir",
                  help="directory with data", metavar="DIR", default="gc-data")
    parser.add_option("--n-type", dest="n_type",
                  help="select n-type", metavar="NTYPE")
    parser.add_option("--sub-category", dest="sub_category",
                  help="select sub-category", metavar="SC")
    parser.add_option("--desc", dest="desc",
                  help="select classifier description", metavar="D")
    parser.add_option("--output-dir", dest="output_dir",
                  help="set output dir", metavar="DIR")
    parser.add_option("--combine", dest="combine", action='store_true',
                  help="combine classifier")
    parser.add_option("--test-classifier", dest="test_classifier", action='store_true',
                  help="test classifier")

    (options, args)  = parser.parse_args()
    return (parser, options, args)

def fn(fid,sc=None,desc=None,n_type=None,opt=None):
    safe_desc = desc
    if safe_desc is not None:
        safe_desc = safe_desc.replace("<","l").replace(">","g")
    if fid=='classifier':
        assert sc is not None
        if desc is None and n_type is None:
            return os.path.join(OUTPUT_DIR, "..", "classifier.%s.json.gz" % sc)
        elif desc is None and n_type is not None:
            return os.path.join(OUTPUT_DIR, "%s/%s/cl_part.json.gz" % (sc,n_type))
        elif desc is not None and n_type is not None:
            return os.path.join(OUTPUT_DIR, "%s/%s/%s/cl_part.json.gz" % (sc,n_type,safe_desc))
        else:
            raise Exception("unknown fid=%s,sc=%s,desc=%s,n_type=%s"%(fid,sc,desc,n_type))
    elif fid in ['doublets-json','doublets-pdb','doublets-params']:
        assert sc is not None
        f = "%s" % sc
        if n_type is not None:
            f += "/" + n_type
        if desc is not None:
            f += "/" + safe_desc
        if opt is not None:
            f += "/doublets-" + opt
            
        if fid=='doublets-json':
            f += ".json.gz"
        elif fid=='doublets-pdb':
            f += ".pdb.gz"
        elif fid=='doublets-params':
            f += ".pickle"
        else:
            raise Exception("unknown fid=%s" % fid)
        return os.path.join(OUTPUT_DIR,f)
    else:
        raise Exception("unknown fid=%s" % fid)

def _reverse_ids(ids):
    res = []
    for id in ids:
        tmp = id.split(":")
        res.append(":".join([tmp[0],tmp[2],tmp[1]]))
    return res

########################

def load_doublets(sc,desc,n_type,options):
    assert len(n_type)==2
    groups = load_json(options.groups)
    rev_n_type = n_type[::-1]
    if desc is not None:
        desc2 = None
        rev_desc2 = None
        if re.match('^[SWH]{2}$',desc):
            rev_desc = desc[::-1]
        else:
            rev_desc = DoubletDescTool.reverse_desc(desc)
            if "cis" in desc:
                desc2 = desc.replace("cis","tran")
                rev_desc2 = rev_desc.replace("cis","tran")
            else:
                desc2 = desc.replace("tran","cis")
                rev_desc2 = rev_desc.replace("tran","cis")
    else:
        rev_desc = None
    all = set(groups['all/all/'+n_type]+_reverse_ids(groups['all/all/'+rev_n_type]))
    ref = set()
    ref2 = set()
    all_ref = set()
    fuzzy = set()
    all_fuzzy = set()
    other = set()
    for k,v in sorted(groups.items()):
        if desc is None:
            if re.match('^classifier/'+sc+'/.*/'+n_type+'$',k):
                ref.update(groups[k])
            if re.match('^classifier/'+sc+'/.*/'+rev_n_type+'$',k):
                ref.update(_reverse_ids(groups[k]))
            if re.match('^fuzzy/'+sc+'/.*/'+n_type+'$',k):
                fuzzy.update(groups[k])
            if re.match('^fuzzy/'+sc+'/.*/'+rev_n_type+'$',k):
                fuzzy.update(_reverse_ids(groups[k]))
        else:
            if re.match('^classifier/'+sc+'/.*/'+n_type+'$',k):
                all_ref.update(groups[k])
            if re.match('^classifier/'+sc+'/.*/'+rev_n_type+'$',k):
                all_ref.update(_reverse_ids(groups[k]))
            if re.match('^fuzzy/'+sc+'/.*/'+n_type+'$',k):
                all_fuzzy.update(groups[k])
            if re.match('^fuzzy/'+sc+'/.*/'+rev_n_type+'$',k):
                all_fuzzy.update(_reverse_ids(groups[k]))
            if re.match('^classifier/'+sc+'/'+desc+'.*/'+n_type+'$',k):
                ref.update(groups[k])
            if re.match('^classifier/'+sc+'/'+rev_desc+'.*/'+rev_n_type+'$',k):
                ref.update(_reverse_ids(groups[k]))

            if desc2 is not None:
                if re.match('^classifier/'+sc+'/'+desc2+'.*/'+n_type+'$',k):
                    ref2.update(groups[k])
                if re.match('^classifier/'+sc+'/'+rev_desc2+'.*/'+rev_n_type+'$',k):
                    ref2.update(_reverse_ids(groups[k]))

            if re.match('^fuzzy/'+sc+'/'+desc+'.*/'+n_type+'$',k):
                fuzzy.update(groups[k])
            if re.match('^fuzzy/'+sc+'/'+rev_desc+'.*/'+rev_n_type+'$',k):
                fuzzy.update(_reverse_ids(groups[k]))
    if desc is None:
        unclassified = all.difference(ref).difference(fuzzy).difference(other)
    else:
        other = all_ref.difference(ref).difference(ref2)
        unclassified = all.difference(ref).difference(ref2).difference(all_fuzzy).difference(other)
    return (ref,fuzzy,other,unclassified)

def compute_doublet_params(ids, params_type, options):
    dp = __import__("doublet-params")

    dd = DoubletsDict(options.data_dir, reduced_atoms=R_ATOMS)
    dd.load_pdb_files(ids, verbose=True)
    res = {}
    
    widgets = ['compute doublet params',' ', Percentage(), ' ', Bar(), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=max(1,len(ids))).start()
    
    for i,d_id in enumerate(ids):
        d = dd.get(d_id)
        if d is None or d[0] is None or d[1] is None:
            print "INVALID doublet! %s" % d_id
            continue
        n_type = dd.get_n_type(d_id).upper()
        if not re.match('^[ACGU]{2}$',n_type):
            print "INVALID doublet! %s, wrong n_type: %s" % (d_id,n_type)
            continue
        d_norm = normalize_points(d, n_type[0])
        d = doublet_params_dict(d, n_type, params_type)
        if d is None:
            print "INVALID doublet! can't compute params %s" % (d_id)
            continue
        c = center_vector(d_norm[1], n_type[1])
        if c is None:
            print "INVALID doublet! missing center %s" % (d_id)
            continue
        d['center'] = c.tolist()
        d['norm2_center'] = c.tolist()
        res[d_id] = d
        # print i,len(ids),d_id,res[d_id]
        pbar.update(i)
    pbar.finish()
    return res

def find_outliers(ref,other,count=100):
    """find indexes of COUNT elements of other closests to ref"""
    if len(ref)==0 or len(other)==0:
        return []
    other_t = KDTree(other)
    distances,indexes = other_t.query(ref,k=max(2,count//10))
    assert len(distances)==len(indexes)
    d = {}
    for i in xrange(len(distances)):
        for j in xrange(len(distances[i])):
            dist,idx = distances[i][j], indexes[i][j]
            if d.has_key(idx):
                d[idx] = min(d[idx], dist)
            else:
                d[idx] = dist
    o = sorted(d.keys(), key=lambda x: d[x])[0:count]
    return set(sorted(d.keys(), key=lambda x: d[x])[0:count])

def remove_outliers(ref,s,count=100):
    if len(s)>1:
        outliers = find_outliers(ref,s,count=count)
        new_s = [x for i,x in enumerate(s) if i not in outliers]
        return new_s
    else:
        return s

def save_points_as_pdb(filename,points,n_type='CG'):
    gp = __import__('gen-pdb')
    s = PDB.Structure.Structure("points")
    for i,p in enumerate(points):
        s.add(gp.points2model(i+1,({},{'C':p}),n_type))
    save_pdb(filename,s)

def save_doublet_points_as_pdb(filename,points,n_type='CG'):
    gp = __import__('gen-pdb')
    s = PDB.Structure.Structure("doublet")
    s.add(gp.points2model(1,points,n_type))
    save_pdb(filename,s)

def save_doublets_as_pdb(filename,doublets,type='bp'):
    gp = __import__('gen-pdb')
    class Opt: pass
    opt = Opt()
    opt.limit = 1000
    opt.type = type
    opt.output = filename
    opt.data_dir = "gc-data"
    opt.add_normalized_base = False
    opt.add_stacking_convex = False
    opt.dont_normalize = False
    opt.quiet = True
    opt.save_json = False
    opt.translate_doublet_ids = False
    opt.only_chain_b = False
    opt.skip_empty = False
    gp.gen_pdb(list(doublets),opt)

def greedy_set_cover(u,sets):
    res = []

    active_sets = set(range(len(sets)))
    for i in active_sets:
        sets[i] = set(sets[i]).intersection(u)

    uncovered = set(u)
    while len(uncovered)>0:
        best_set = None
        best_len = 0
        for i in active_sets:
            if len(sets[i])>best_len:
                best_set = i
                best_len = len(sets[i])
        if best_set is None:
            break
        res.append(best_set)
        active_sets.remove(best_set)
        uncovered.difference_update(sets[best_set])
        for i in active_sets:
            sets[i].difference_update(sets[best_set])
    return (res, uncovered)

# @cache("compute_rmsd_distances")
def compute_rmsd_distances(dd,ref,other):
    dist = []
    widgets = ['compute rmsd distances',' ', Percentage(), ' ', Bar(), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=max(1,len(ref))).start()

    other_norm = [dd.get_normalized(x) for x in other]
    for i,d_id in enumerate(ref):
        pbar.update(i)
        n_type = dd.get_n_type(d_id)
        p = dd.get_normalized(d_id)
        dist.append([bp_distance(x,p,n_type,already_normalized=True) for x in other_norm])
    pbar.finish()
 
    return dist


def rmsd_classifier(ref,other,uncl,options,outliers_count=100,cl_count=10,dist_factor=1.0):    
    org_ref = [x for x in ref]
    if len(ref)>1000:
        random.shuffle(ref)
        ref = ref[0:1000]

    ids = set(ref+other+uncl)
    dd = DoubletsDict(options.data_dir, reduced_atoms=R_ATOMS)
    dd.load_pdb_files(ids, verbose=True)

    dist_other = compute_rmsd_distances(dd,ref,other)
    dist_uncl = compute_rmsd_distances(dd,ref,uncl)
    dist_ref = compute_rmsd_distances(dd,ref,ref)

    print "rmsd classifier, before remove outliers", len(other), len(uncl)

    tmp_uncl = [100.0]*len(uncl)
    tmp_other = [100.0]*len(other)
    for i,d_id in enumerate(ref):
        for j,dd_id in enumerate(other):
            tmp_other[j] = min(tmp_other[j], dist_other[i][j])
        for j,dd_id in enumerate(uncl):
            tmp_uncl[j] = min(tmp_uncl[j], dist_uncl[i][j])
    removed_other = 0
    removed_uncl = 0
    for j in sorted(range(len(other)),key=lambda x:tmp_other[x])[0:min(len(other)//10, outliers_count)]:
        removed_other += 1
        for i,d_id in enumerate(ref):
            dist_other[i][j]=100.0
    for j in sorted(range(len(uncl)),key=lambda x:tmp_uncl[x])[0:min(len(uncl)//10, outliers_count)]:
        removed_uncl += 1
        for i,d_id in enumerate(ref):
            dist_uncl[i][j]=100.0

    print "rmsd classifier, before remove outliers", len(other)-removed_other, len(uncl)-removed_uncl

    d = []
    for i,d_id in enumerate(ref):
        d.append(min(min(dist_other[i]+[8.0])/2, min(dist_uncl[i]+[8.0])))
        # print i,d_id,"%.3f"%d[i]

    if len(ref)<=cl_count:
        selected_ids = range(len(ref))
        uncovered = []
    else:
        u = xrange(len(ref))
        sets = [[j for j in xrange(len(ref)) if dist_ref[i][j]<=d[i]] for i in xrange(len(ref))]
        res, uncovered = greedy_set_cover(u,sets)
        selected_ids = res[0:cl_count]
        if len(selected_ids)<cl_count:
            unselected = [i for i in xrange(len(ref)) if i not in selected_ids]
            random.shuffle(unselected)
            selected_ids += unselected[0:(cl_count-len(selected_ids))]
            
    result = [(d[i]*dist_factor,ref[i]) for i in selected_ids]
    print "rmsd_classifier", result, len(uncovered)
    return result

def compute_distance_matrix(ref,options):
    org_ref = [x for x in ref]
    if len(ref)>100:
        random.shuffle(ref)
        ref = ref[0:100]

    ids = set(ref)
    dd = DoubletsDict(options.data_dir, reduced_atoms='*')
    dd.load_pdb_files(ids, verbose=True)
    
    print "computing distance matrix"

    dm = {}
    for d_id in ref:
        (p1,p2) = dd.get(d_id)
        for k1,a1 in p1.items():
            p1[k1] = np.array(a1,'f')
        for k2,a2 in p2.items():
            p2[k2] = np.array(a2,'f')
        for k1,a1 in p1.items():
            if not dm.has_key(k1):
                dm[k1] = {}
            for k2,a2 in p2.items():
                d = vector_length(a1-a2)
                d = float(d)
                if not dm[k1].has_key(k2):
                    dm[k1][k2] = {"min_d":d,"max_d":d,"sum_d":0.0,"c":0}
                dm[k1][k2]['sum_d']+=d
                dm[k1][k2]['c']+=1
    for k1,d_row in dm.items():
        for k2 in d_row.keys():
            dm[k1][k2]['avg_d'] = dm[k1][k2]['sum_d']/dm[k1][k2]['c']
            del dm[k1][k2]['sum_d']
    return dm


def point_classifier(ref,other,uncl,outliers_count=100,cl_count=500,cl_groups=10,dist_factor=0.5):
    other = remove_outliers(ref,other,count=outliers_count)
    uncl = remove_outliers(ref,uncl,count=outliers_count)
    
    if len(other)==0:
        other = np.array([(0,0,0)],'f')
    if len(uncl)==0:
        uncl = np.array([(0,0,0)],'f')

    r,_d = kmeans(ref,cl_count)

    other_t = KDTree(other)
    uncl_t = KDTree(uncl)

    distances_r = []
    for i,p in enumerate(r):
        distances_o,_index_o = other_t.query(p,k=1)
        distances_u,_index_u = uncl_t.query(p,k=1)
        d = min(distances_o, distances_u)
        # print i,d
        distances_r.append(d)
    r_seq = sorted(range(len(r)), key=lambda x: distances_r[x], reverse=True)
    
    classifier = []
    
    for group_num,r_group in enumerate(np.array_split(r_seq, min(cl_groups, len(r_seq)))):
        min_d = min([distances_r[x] for x in r_group])
        max_d = max([distances_r[x] for x in r_group])
        group_dist = min_d + dist_factor * (max_d-min_d)
        group_points = [r[x].tolist() for x in r_group]
        # print len(r_seq), min_d, max_d, group_dist
        # save_points_as_pdb("/tmp/aaa_%s_dist_%.2f.pdb"%(group_num,group_dist),[r[x] for x in r_group])
        classifier.append((group_dist,group_points))

    return classifier

def points_classifier_check(cl, points):
    res = []
    cl_trees = [(cl_d,KDTree(cl_p)) for cl_d,cl_p in cl]
    for i,p in enumerate(points):
        ok = False
        for dist,t in cl_trees:
            d,j = t.query(p,k=1)
            if d<dist:
                # print "OK i=%d d=%.3f" % (i,d)
                ok = True
                break
        if ok:
            res.append(i)
    return res
    

def filter_by_params(doublets, limits, all_params_dict):
    res = []
    for d_id in doublets:
        if not all_params_dict.has_key(d_id):
            continue
        p = all_params_dict[d_id]
        ok = True
        for k,limit in limits.items():
            if limit.has_key('min'):
                if not (p[k]>=limit['min'] and p[k]<=limit['max']):
                    ok = False
                    break
            elif limit.has_key('min1'):
                if not ( (p[k]>=limit['min1'] and p[k]<=limit['max1']) or (p[k]>=limit['min2'] and p[k]<=limit['max2']) ):
                    ok = False
                    break
        if ok:
            res.append(d_id)
    return res


def filter_by_distance(doublets, dist_limits, all_params_dict):
    res = []
    d_trees = [(d,KDTree(d_points)) for d,d_points in dist_limits]
    for i,d_id in enumerate(doublets):
        if not all_params_dict.has_key(d_id):
            continue
        p = all_params_dict[d_id]
        if not p.has_key('center'):
            continue
        ok = False
        for dist,t in d_trees:
            d,j = t.query(p['center'],k=1)
            if d<dist:
                # print "OK i=%d d=%.3f" % (i,d)
                ok = True
                break
        if ok:
            res.append(d_id)
    return res

def filter_points_by_distance(points, dist_limits):
    res = []
    d_trees = [(d,KDTree(d_points)) for d,d_points in dist_limits]
    for i,p in enumerate(points):
        ok = False
        for dist,t in d_trees:
            d,j = t.query(p,k=1)
            if d<dist:
                ok = True
                break
        if ok:
            res.append(p)
    return res


def filter_by_rmsd_distance(doublets, rmsd_limits, dd):
    res = []
    r = [(d,dd.get_normalized(x)) for d,x in rmsd_limits]
    for i,d_id in enumerate(doublets):
        n_type = dd.get_n_type(d_id)
        p = dd.get_normalized(d_id)
        ok = False
        for dist,ref_d in r:
            d = bp_distance(p,ref_d,n_type,already_normalized=True)
            if d<dist:
                ok = True
                break
        if ok:
            res.append(d_id)
    return res

def categorize_by_rmsd_distance(doublets, rmsd_limits, dd):
    res = [[] for i in xrange(len(rmsd_limits))]
    r = [(d,dd.get_normalized(x)) for d,x in rmsd_limits]
    
    widgets = ['categorize by rmsd',' ', Percentage(), ' ', Bar(), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=max(1,len(doublets))).start()
    
    for i,d_id in enumerate(doublets):
        n_type = dd.get_n_type(d_id)
        p = dd.get_normalized(d_id)
        ok = False
        bin_num = None
        min_d = None
        for j,(dist,ref_d) in enumerate(r):
            d = bp_distance(p,ref_d,n_type,already_normalized=True)
            if d<dist:
                if min_d is None or d<min_d:
                    bin_num = j
                    min_d = d
        if bin_num is not None:
            res[bin_num].append(d_id)
        pbar.update(i)
    pbar.finish()
    return res

def compute_bp_classifier(options):
    (ref,fuzzy,_other,unclassified) = load_doublets('bp',None,options.n_type,options)

    print "loading doublets, ref=%d fuzzy=%d unclassified=%d" % (len(ref), len(fuzzy), len(unclassified))

    all_ids = ref.union(fuzzy).union(unclassified)
    all_params = compute_doublet_params(all_ids, 'bp', options)
    pickle_save(fn("doublets-params",sc="bp",n_type=options.n_type),all_params)
        
    for k in 'nn_ang','nn_ang_norm','dist','min_dist','n1cc_ang','n2cc_ang','n12cc_ang':
        values = [all_params[d_id][k] for d_id in ref]
        min_v = min(values)
        max_v = max(values)
        values.sort()
        if len(values)>100:
            min99_v = values[len(values)//100]
            max99_v = values[len(values)-len(values)//100]
        elif len(values)>30:
            min99_v = values[len(values)//30]
            max99_v = values[len(values)-len(values)//30]
        else:
            min99_v = min_v
            max99_v = max_v
        print "k=%s min=%.3f max=%.3f min99%%=%.3f max99%%=%.3f" % (k, min_v, max_v, min99_v, max99_v)

    limits = {'nn_ang_norm': {'min':0.0,'max':65.0}, 
              'dist': {'min':4.0,'max':8.5}, 
              'min_dist': {'min':0.1,'max':3.2}, 
              'n1cc_ang': {'min': 50.0,'max':140.0}, 
              'n2cc_ang': {'min': 50.0,'max':140.0}
              }
    
    ref = filter_by_params(ref, limits, all_params) 
    fuzzy = filter_by_params(fuzzy, limits, all_params) 
    unclassified = filter_by_params(unclassified, limits, all_params) 

    print "after filtering by params, ref=%d fuzzy=%d unclassified=%d" % (len(ref), len(fuzzy), len(unclassified))

    if False:
        d_limits = point_classifier(
                array([all_params[d_id]['center'] for d_id in ref],'f'),
                array([(0,0,0)],'f'),
                array([all_params[d_id]['center'] for d_id in unclassified],'f'),
                outliers_count=4000,
                cl_count=800,
                cl_groups=5,
                dist_factor=1.4
        )
        
        ref = filter_by_distance(ref, d_limits, all_params) 
        fuzzy = filter_by_distance(fuzzy, d_limits, all_params) 
        unclassified = filter_by_distance(unclassified, d_limits, all_params) 
    
        print "after filtering by distance, ref=%d fuzzy=%d unclassified=%d" % (len(ref), len(fuzzy), len(unclassified))
    
    cl3 = ["RETURN","?bp"]
    # cl2 = ["IF",[["OR"]+[["LE",["DISTANCE_TO_POINTS","$norm2_center",p],d] for d,p in d_limits]],cl3]
    cl1 = ["IF",[["IN_RANGE","$d_param_%s"%k,v['min'],v['max']] for k,v in limits.items()],[cl3]]
    cl0 = ["IF",[["EQ","$n_type",options.n_type]],[cl1]]
    save_json(fn("classifier",sc="bp",n_type=options.n_type), [cl0], indent=2)
    
    for suffix,points,doublets in (
        ["ref",[all_params[d_id]['center'] for d_id in ref],ref],
        ["fuzzy",[all_params[d_id]['center'] for d_id in fuzzy],fuzzy],
        ["unclassified",[all_params[d_id]['center'] for d_id in unclassified],unclassified]
    ):
        save_points_as_pdb(fn("doublets-pdb",sc="bp",n_type=options.n_type,opt=suffix), points, options.n_type)
        save_json(fn("doublets-json",sc="bp",n_type=options.n_type,opt=suffix), doublets)

def compute_stacking_classifier(options):
    pass

def join_classifiers(cl1, cl2):
    def walk(x):
        if isinstance(x,list):
            if len(x)==1 and isinstance(x[0],list) and len(x[0])==2 and x[0][0]=='RETURN':
                return cl2
            else:
                new_res = []
                for el in x:
                    new_res.append(walk(el))
                return new_res
        else:
            return x
    return walk(cl1)
    return cl1
    
def classifier_bp_cis_tran(desc,n_type,cl_cis,cl_tran):
    cis_value = expected_strand_orient(desc[0:2]+"_cis")
    tran_value = expected_strand_orient(desc[0:2]+"_tran")
    return [
                ["IF",[["EQ","$d_param_strand_orient",cis_value]],cl_cis],
                ["IF",[["EQ","$d_param_strand_orient",tran_value]],cl_tran],
            ]

def limits_condition(limits,ok_code):
    cond = []
    for key,limit in limits.items():
        if limit.has_key('min1'):
            cond.append(["OR",
                        ["IN_RANGE","$d_param_%s"%key,limit['min1'],limit['max1']],
                        ["IN_RANGE","$d_param_%s"%key,limit['min2'],limit['max2']],
                   ])
        else:
            cond.append(["IN_RANGE","$d_param_%s"%key,limit['min'],limit['max']])
    return ["IF",cond,ok_code]

def compute_bp_spec_classifier(options):
    (ref,fuzzy,other,unclassified) = load_doublets('bp',options.desc,options.n_type,options)
    status_msg("loading doublets, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified)))

    ref = list(set(ref).intersection(set(load_json(fn("doublets-json",sc="bp",n_type=options.n_type,opt="ref")))))
    if len(ref)==0:
        status_msg("empty reference set!")
        save_json(fn("classifier",sc="bp",desc=options.desc,n_type=options.n_type), [], indent=2)
        return
    other = list(set(other).intersection(set(load_json(fn("doublets-json",sc="bp",n_type=options.n_type,opt="ref")))))
    fuzzy = list(set(fuzzy).intersection(set(load_json(fn("doublets-json",sc="bp",n_type=options.n_type,opt="fuzzy")))))
    unclassified = list(set(unclassified).intersection(set(load_json(fn("doublets-json",sc="bp",n_type=options.n_type,opt="unclassified")))))
    
    ref.sort()
    other.sort()
    fuzzy.sort()
    unclassified.sort()

    status_msg("loading doublets, after filtering, ref=%d fuzzy=%d unclassified=%d" % (len(ref), len(fuzzy), len(unclassified)))
    
    all_params = pickle_load(fn("doublets-params",sc="bp",n_type=options.n_type))
    assert all_params is not None
    
    limits = None
    strict_limits = None
    if "cis" in options.desc or "tran" in options.desc:
        ee = expected_strand_orient(options.desc)
        limits = {"strand_orient":{"min":ee,"max":ee}}
        strict_limits = {}
        eps = 0.001
        
        import bp_limits
        
        for key in ['nn_ang_norm','min_dist','rot_ang','dist']:
            pr = "99.00"
            if (options.desc,options.n_type) in [
                    ('WW_cis','CG'), ('WW_cis','GC'), ('WW_cis','AU'), ('WW_cis','UA'), ('WW_cis','UG'), ('WW_cis','GU'),
                    ('SH_tran','GA'), ('HS_tran','AG'), ('WH_tran','UA'), ('HW_tran','AU'),
            ]:
                pr = "99.90"
            elif key=='dist':
                pr = "99.90"
            l = bp_limits.limits.get('%s/%s/%s/%s' % (options.desc, options.n_type, key, pr))
            if l is not None:
                if key in ['nn_ang_norm']:
                    l['min'] = 0.0
                for k in ['min','min1','min2']:
                    if l.has_key(k):
                        l[k] = l[k]-eps
                for k in ['max','max1','max2']:
                    if l.has_key(k):
                        l[k] = l[k]+eps
                strict_limits[key] = l
                print "param=%s limits=%s" % (key, strict_limits[key])
        
        ref = filter_by_params(ref, limits, all_params) 
        other = filter_by_params(other, limits, all_params) 
        fuzzy = filter_by_params(fuzzy, limits, all_params) 
        unclassified = filter_by_params(unclassified, limits, all_params) 
        status_msg("after orient filtering, ref=%d other=%d fuzzy=%d unclassified=%d" % (len(ref), len(other), len(fuzzy), len(unclassified)))

        ref = filter_by_params(ref, strict_limits, all_params) 
        other = filter_by_params(other, strict_limits, all_params) 
        fuzzy = filter_by_params(fuzzy, strict_limits, all_params) 
        unclassified = filter_by_params(unclassified, strict_limits, all_params) 

        status_msg("after strict limits filtering, ref=%d other=%d fuzzy=%d unclassified=%d" % (len(ref), len(other), len(fuzzy), len(unclassified)))

    
    points = [all_params[d_id]['center'] for d_id in ref]
    save_points_as_pdb(fn("doublets-pdb",sc="bp",desc=options.desc,n_type=options.n_type,opt="ref-centers"), points, options.n_type)
    save_json(fn("doublets-json",sc="bp",desc=options.desc,n_type=options.n_type,opt="ref-centers"), ref)
    points = [all_params[d_id]['center'] for d_id in other]
    save_points_as_pdb(fn("doublets-pdb",sc="bp",desc=options.desc,n_type=options.n_type,opt="other-centers"), points, options.n_type)
    save_json(fn("doublets-json",sc="bp",desc=options.desc,n_type=options.n_type,opt="other-centers"), other)
    points = [all_params[d_id]['center'] for d_id in unclassified]
    save_points_as_pdb(fn("doublets-pdb",sc="bp",desc=options.desc,n_type=options.n_type,opt="uncl-centers"), points, options.n_type)
    save_json(fn("doublets-json",sc="bp",desc=options.desc,n_type=options.n_type,opt="uncl-centers"), other)

    if len(ref)==0:
        status_msg("empty reference set!")
        save_json(fn("classifier",sc="bp",desc=options.desc,n_type=options.n_type), [], indent=2)
        return
    
    cl_count=300
    dist_factor=1.2
    cl_groups=5
    outliers_count=1000
    if (options.desc,options.n_type) in (
            ('SH_tran','GA'),('HS_tran','AG'),
            ('SS_tran','AG'),('SS_tran','GA'),
            ('SS_cis','GA'),('SS_cis','AG'),
            ('WW_cis','CG'),('WW_cis','GC'),
            ('HW_tran','UA'),('WH_tran','AU')):
        cl_count=800
        dist_factor=1.4
        outliers_count=2000
        cl_groups=8
    d_limits = point_classifier(
        array([all_params[d_id]['center'] for d_id in ref],'f'),
        # array([all_params[d_id]['center'] for d_id in other],'f'),
        array([(0,0,0)],'f'),
        array([all_params[d_id]['center'] for d_id in unclassified],'f'),
        outliers_count=outliers_count,
        cl_count=cl_count,
        cl_groups=cl_groups,
        dist_factor=dist_factor
    )

    ref = filter_by_distance(ref, d_limits, all_params) 
    other = filter_by_distance(other, d_limits, all_params) 
    fuzzy = filter_by_distance(fuzzy, d_limits, all_params) 
    unclassified = filter_by_distance(unclassified, d_limits, all_params) 

    status_msg("after filter by distance, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified)))

    save_doublets_as_pdb(fn("doublets-pdb",sc="bp",desc=options.desc,n_type=options.n_type,opt="ref"), ref, "bp")
    save_json(fn("doublets-json",sc="bp",desc=options.desc,n_type=options.n_type,opt="ref"), ref)
    
    #points = [all_params[d_id]['center'] for d_id in other]
    save_doublets_as_pdb(fn("doublets-pdb",sc="bp",desc=options.desc,n_type=options.n_type,opt="other"), other, "bp")
    save_json(fn("doublets-json",sc="bp",desc=options.desc,n_type=options.n_type,opt="other"), other)
    
    dist_factor = 1.2
    outliers_count = 16
    cl_count = 16
    if (options.desc,options.n_type) in (
            # ('SH_tran','GA'),('HS_tran','AG'),
            ('SS_tran','GA'),('SS_tran','AG'),
            ('SS_cis','GA'),('SS_cis','AG'),
            ('HW_tran','UA'),('WH_tran','AU'),
            ('HW_tran','GU'),('WH_tran','UG'),
            ('HH_cis','CG'),('HH_cis','GC'),
            # dodane 2013-01-28
            ('HW_cis','AC'),('WH_cis','CA'),
            ('WS_tran','AC'),('SW_tran','AC'),
            ('HW_tran','CU'),('WH_tran','UC'),
            ('WH_tran','CU'),('HW_tran','UC'),
            ('SW_cis','CU'),('SW_cis','UC'),
            ('WS_cis','UU'),('SW_cis','UU'),
            ('HS_cis','UU'),('SH_cis','UU'),
            ('HS_cis','AG'),('SH_cis','GA'),
        ):
        dist_factor = 1.5
        outliers_count = 100
    elif (options.desc,options.n_type) in (
            ('SH_tran','GA'),('HS_tran','AG'),
        ):
        dist_factor = 1.3
        outliers_count = 50

    # [2013-02-25] removed: options.n_type in ('AA','CC','GG','UU')
    if (options.desc,options.n_type) in (
            ('WW_cis','CG'),('WW_cis','GC'),
            ('WW_cis','AU'),('WW_cis','UA'),
            ('SS_cis','AG'),('SS_cis','GA'),
            ('SS_tran','AG'),('SS_tran','GA'),
        ):
        cl_count=32
        
    status_msg("computing rmsd classifier")

    rmsd_limits = rmsd_classifier(
            ref, other, unclassified,
            options=options,
            outliers_count=outliers_count,
            cl_count=cl_count,
            dist_factor=dist_factor
    )

    status_msg("computing rmsd classifier - done")
    
    ids = set(ref+other+fuzzy+unclassified)
    dd = DoubletsDict(options.data_dir,reduced_atoms=R_ATOMS)
    dd.load_pdb_files(ids, verbose=True)
    
    ref = filter_by_rmsd_distance(ref, rmsd_limits, dd)
    other = filter_by_rmsd_distance(other, rmsd_limits, dd)
    fuzzy = filter_by_rmsd_distance(fuzzy, rmsd_limits, dd)
    unclassified = filter_by_rmsd_distance(unclassified, rmsd_limits, dd)

    status_msg("after filtering by rmsd, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified)))

    ref_cat = categorize_by_rmsd_distance(ref, rmsd_limits, dd)
    assert len(rmsd_limits)==len(ref_cat)
    # organize by size of groups
    idx = sorted(range(len(ref_cat)), key=lambda x: len(ref_cat[x]), reverse=True)
    rmsd_limits = [rmsd_limits[i] for i in idx]
    ref_cat = [ref_cat[i] for i in idx]
    cl_doublets = [d_id for (l,d_id) in rmsd_limits]
    save_doublets_as_pdb(fn("doublets-pdb",sc="bp",desc=options.desc,n_type=options.n_type,opt="cl"), cl_doublets, "full")
    save_json(fn("doublets-json",sc="bp",desc=options.desc,n_type=options.n_type,opt="cl"), cl_doublets)

    dm_cl = None
    for j,d in enumerate(ref_cat,start=1):
        status_msg("saving reference doublets-%d count=%d" % (j,len(d)))
        save_doublets_as_pdb(fn("doublets-pdb",sc="bp",desc=options.desc,n_type=options.n_type,opt="ref-%02d"%j), d[0:1000], "bp")
        save_json(fn("doublets-json",sc="bp",desc=options.desc,n_type=options.n_type,opt="ref-%02d"%j), d)
        if j==1 and len(d)>500: 
            dm = compute_distance_matrix(d,options)
            dm_cl = []
            dm_cl.append(["COMPARE_DISTANCES", dm, d[0]])
            dm_cl.append(["RETURN_IF", "$dist_score", 0.0001, 0.01, "%s (dm)" %(options.desc)])


    status_msg("saving classifier")
    cl1 = []
    for i,(d,ref) in enumerate(rmsd_limits,start=1):
        r1,r2 = dd.get_normalized(ref)
        for k,v in r1.items():
            r1[k] = v.tolist()
        for k,v in r2.items():
            r2[k] = v.tolist()
        cl1.append(["COMPUTE_BP_DISTANCE",(r1,r2),ref])
        if "cis" in options.desc or "tran" in options.desc:
            if strict_limits is not None:
                cl1 += [limits_condition(strict_limits, [["RETURN_IF","$bp_distance",d,d,"%s (%d)"%(options.desc,i)]])]
                cl1 += [["RETURN_IF","$bp_distance",0.0001,2*d,"%s (f%d)"%(options.desc,i)]]
            else:
                cl1 += [["RETURN_IF","$bp_distance",d,2*d,"%s (%d)"%(options.desc,i)]]
        else:
            # TODO: usunac ten kod, juz nie uzywamy tego typu klasyfikatora!
            ret_cis = [["RETURN_IF","$bp_distance",d,2*d,"%s_cis (%d)"%(options.desc,i)]]
            ret_tran = [["RETURN_IF","$bp_distance",d,2*d,"%s_tran (%d)"%(options.desc,i)]]
            cl1 += classifier_bp_cis_tran(options.desc,options.n_type,ret_cis,ret_tran)
    cl2 = ["IF",[["OR"]+[["LE",["DISTANCE_TO_POINTS","$norm2_center",p],d] for d,p in d_limits]],cl1]

    if limits is not None:
        cl2 = ["IF",[["IN_RANGE","$d_param_%s"%k,v['min'],v['max']] for k,v in limits.items()],[cl2]]

    cl_full = [cl2]
    if dm_cl is not None:
        cl_full += dm_cl

    save_json(fn("classifier",sc="bp",desc=options.desc,n_type=options.n_type), cl_full, indent=2)
    status_msg("saving classifier - done")

def compute_stacking_spec_classifier(options):
    (ref,fuzzy,other,unclassified) = load_doublets('stacking',options.desc,options.n_type,options)
    print "loading doublets, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified))

    all_ids = ref.union(fuzzy).union(unclassified)
    all_params = compute_doublet_params(all_ids, 'stacking', options)
    pickle_save(fn("doublets-params",sc="stacking",n_type=options.n_type),all_params)
        
    for k in 'nn_ang_norm','dist','n12cc_ang','stack_orient':
        p = sorted([all_params[d_id][k] for d_id in ref])
        min_p = None
        max_p = None
        min_p99 = None
        max_p99 = None
        if len(p)>0:
            min_p = min(p)
            max_p = max(p)
        if len(p)>100:
            min_p99 = p[len(p)//100]
            max_p99 = p[len(p)-len(p)//100]
        elif len(p)>30:
            min_p99 = p[len(p)//30]
            max_p99 = p[len(p)-len(p)//30]
        else:
            min_p99 = min(p)
            max_p99 = max(p)
        if len(p)>0:
            print "param: %s min=%.2f max=%.2f min(99%%)=%.2f max(99%%)=%.2f" % (k,min_p,max_p,min_p99,max_p99)

    common_limits = {'stack_orient': {'min': expected_stack_orient(options.desc), 'max': expected_stack_orient(options.desc)},
                    'nn_ang_norm': { 'min':0.0,'max':45.0}, 
                    'dist': {'min':0.0,'max':6.5}, 
                    'n12cc_ang': {'min': 0.0,'max':60.0} }
                    
    fr3d_limits = {
                    'stack_overlap': {'min': 1.0,'max':100.0},
                    'min_dist': {'min': 1.0,'max':4.0},
                    'stack_norm': {'min': 0.6,'max':999.0},
                  }
    
    ref = filter_by_params(ref,common_limits,all_params)
    other = filter_by_params(other,common_limits,all_params)
    fuzzy = filter_by_params(fuzzy,common_limits,all_params)
    unclassified = filter_by_params(unclassified,common_limits,all_params)

    print "after filtering by common params, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified))


    cl = []
    for name,limits in [
        ('mo',{'nn_ang_norm': {'min':0.0,'max':30.0}, 
              'dist': {'min':0.0,'max':5.5}, 
              'n12cc_ang': {'min': 0.0,'max':40.0} }),
        ('relaxed-dist', {'nn_ang_norm': {'min':0.0,'max':30.0}, 
              'dist': {'min':5.5,'max':6.0}, 
              'n12cc_ang': {'min': 0.0,'max':40.0} }),
        ('relaxed-nn', {'nn_ang_norm': {'min':30.0,'max':38.0}, 
              'dist': {'min':0.0,'max':6.0}, 
              'n12cc_ang': {'min': 0.0,'max':40.0} }),
        ('relaxed-n12cc', {'nn_ang_norm': {'min':0.0,'max':30.0}, 
              'dist': {'min':0.0,'max':6.0}, 
              'n12cc_ang': {'min': 40.0,'max':50.0} }),
        ('more-relaxed-n12cc', {'nn_ang_norm': {'min':0.0,'max':30.0}, 
              'dist': {'min':0.0,'max':6.0}, 
              'n12cc_ang': {'min': 50.0,'max':55.0} }),
        ('fuzzy', {'nn_ang_norm': { 'min':0.0,'max':45.0}, 
              'dist': {'min':0.0,'max':6.5}, 
              'n12cc_ang': {'min': 0.0,'max':60.0} }), 
        ]:
        if name in ['relaxed-nn','relaxed-n12cc','more-relaxed-n12cc']:
            for k,v in fr3d_limits.items():
                limits[k]=v
        
        cur_ref = filter_by_params(ref,limits,all_params)
        cur_other = filter_by_params(other,limits,all_params)
        cur_fuzzy = filter_by_params(fuzzy,limits,all_params)
        cur_unclassified = filter_by_params(unclassified,limits,all_params)
        print "%s: after filtering params, ref=%d fuzzy=%d other=%d unclassified=%d" % (name, len(cur_ref), len(cur_fuzzy), len(cur_other), len(cur_unclassified))

        points = [all_params[d_id]['center'] for d_id in cur_ref]
        save_points_as_pdb(fn("doublets-pdb",sc="stacking",desc=options.desc,n_type=options.n_type,opt="ref-centers-%s"%name), points, options.n_type)
        save_json(fn("doublets-json",sc="stacking",desc=options.desc,n_type=options.n_type,opt="ref-centers"), cur_ref)
        points = [all_params[d_id]['center'] for d_id in cur_other]
        save_points_as_pdb(fn("doublets-pdb",sc="stacking",desc=options.desc,n_type=options.n_type,opt="other-centers-%s"%name), points, options.n_type)
        save_json(fn("doublets-json",sc="stacking",desc=options.desc,n_type=options.n_type,opt="other-centers"), cur_other)
        points = [all_params[d_id]['center'] for d_id in cur_unclassified]
        save_points_as_pdb(fn("doublets-pdb",sc="stacking",desc=options.desc,n_type=options.n_type,opt="uncl-centers-%s"%name), points, options.n_type)
        save_json(fn("doublets-json",sc="stacking",desc=options.desc,n_type=options.n_type,opt="uncl-centers"), cur_other)

        d_limits = None
        if len(cur_ref)>0:
            dist_factor = 1.0
            if name in ['fuzzy']:
                dist_factor = 1.2
            elif name in ['mo']:
                dist_factor = 3.0
            d_limits = point_classifier(
                array([all_params[d_id]['center'] for d_id in cur_ref],'f'),
                array([(0,0,0)],'f'),
                array([all_params[d_id]['center'] for d_id in unclassified],'f'),
                outliers_count=8,
                cl_count=1600,
                cl_groups=16,
                dist_factor=dist_factor
            )
        if d_limits is not None:
            cur_ref = filter_by_distance(cur_ref, d_limits, all_params) 
            cur_other = filter_by_distance(cur_other, d_limits, all_params) 
            cur_fuzzy = filter_by_distance(cur_fuzzy, d_limits, all_params) 
            cur_unclassified = filter_by_distance(cur_unclassified, d_limits, all_params) 
            print "%s: after filtering by distance, ref=%d fuzzy=%d other=%d unclassified=%d" % (name, len(cur_ref), len(cur_fuzzy), len(cur_other), len(cur_unclassified))

            
        save_doublets_as_pdb(fn("doublets-pdb",sc="stacking",desc=options.desc,n_type=options.n_type,opt="ref-%s"%name), cur_ref, "bp")
        save_json(fn("doublets-json",sc="stacking",desc=options.desc,n_type=options.n_type,opt="ref-%s"%name), cur_ref)
        
        if len(cur_ref)>0 and name in ['mo']:
            cl_doublets = [x for x in cur_ref]
            random.shuffle(cl_doublets)
            cl_doublets = cl_doublets[0:8]
            save_doublets_as_pdb(fn("doublets-pdb",sc="stacking",desc=options.desc,n_type=options.n_type,opt="cl"), cl_doublets, "full")
            save_json(fn("doublets-json",sc="stacking",desc=options.desc,n_type=options.n_type,opt="cl"), cl_doublets)
        
        if len(cur_ref)>0 or name in ['mo']:
            cl_part = [["RETURN","%s (%s)"%(options.desc,name)]]
            if d_limits is not None:
                cl_part = []
                for num,(d,p) in enumerate(d_limits,start=1):
                    if name in ['fuzzy']:
                        dd = 0.001
                    else:
                        dd = d
                    dd2 = 2*d
                    cl_part += [["COMPUTE_DISTANCE_TO_POINTS","$norm2_center",p]]
                    cl_part += [["RETURN_IF","$distance_to_points",dd,dd2,"%s (%s-%d)"%(options.desc,name,num)]]
            if limits is not None:
                cl_part = [["IF",[["IN_RANGE","$d_param_%s"%k,v['min'],v['max']] for k,v in limits.items()],cl_part]]
        
            cl += cl_part
        else:
            print "%s: ignoring, empty reference set" % (name)
        
    
    cl = ["IF",[["EQ","$n_type",options.n_type]]+[["IN_RANGE","$d_param_%s"%k,v['min'],v['max']] for k,v in common_limits.items()],cl]

    save_json(fn("classifier",sc="stacking",desc=options.desc,n_type=options.n_type), [cl], indent=2)

def combine_classifier(options):
    assert options.sub_category is not None
    if options.sub_category=='bp':
        DESCS = []
        DESCS += [x1+x2+"_cis" for x1,x2 in itertools.product(['W','S','H'],repeat=2)]
        DESCS += [x1+x2+"_tran" for x1,x2 in itertools.product(['W','S','H'],repeat=2)]
    elif options.sub_category=='stacking':
        DESCS = ['<<','<>','><','>>']
    elif options.sub_category=='base-phosphate':
        DESCS = ["H_0BPh", "SW_2BPh", "S_1BPh", "W_345BPh", "W_6BPh", "H_789BPh"]
    elif options.sub_category=='base-ribose':
        DESCS = ["H_0BR", "SW_2BR", "S_1BR", "W_345BR", "W_6BR", "H_789BR"]
    elif options.sub_category=='other':
        DESCS = ["diagonal-c","diagonal-nc-ww","long-stacking-c"]
    elif options.sub_category=='other2':
        DESCS = ['base-ribose-stacking']
    elif options.sub_category=='other3':
        DESCS = CL_CAT[options.sub_category]
    else:
        raise Exception("unsupported sub_category=%s" % options.sub_category)
    res = []
    for n_type in N_TYPES:
        if options.sub_category in ['bp']:
            cl1 = load_json(fn("classifier",sc=options.sub_category,n_type=n_type))
        cl2_sum = []
        for desc in DESCS:
            f = fn("classifier",sc=options.sub_category,desc=desc,n_type=n_type)
            if os.path.isfile(f):
                cl2 = load_json(f)
                if not isinstance(cl2,list):
                    print "corrupted file: %s" % f
                    continue
                if len(cl2)>0:
                    print "OK: %s" % f
                    cl2_sum += cl2
                else:
                    print "empty classifier: %s" % f
            else:
                print "missing: %s" % f
        if options.sub_category in ['bp']:
            res += join_classifiers(cl1,cl2_sum)
        else:
            res += cl2_sum
    save_json(fn("classifier",sc=options.sub_category),res,indent=2)

def bph_and_br_fr3d_classifier(options):
    """
    %  1  A C2-H2  interacts with oxygen of phosphate, called 2BPh
    %  2  A N6-1H6 interacts with oxygen of phosphate, called 6BPh
    %  3  A N6-2H6 interacts with oxygen of phosphate, called 7BPh
    %  4  A C8-H8  interacts with oxygen of phosphate, called 0BPh
    %
    %  5  C N4-2H4 interacts with oxygen of phosphate, called 6BPh
    %  6  C N4-1H4 interacts with oxygen of phosphate, called 7BPh
    %  7  C N4-1H4 and C5-H5 interact with 2 oxygens of phosphate, called 8BPh
    % 18  C N4-1H4 and C5-H5 interact with just one oxygen, called 7BPh
    %  8  C C5-H5  interacts with oxygen of phosphate, called 9BPh
    %  9  C C6-H6  interacts with oxygen of phosphate, called 0BPh
    %
    % 10  G N2-1H2 interacts with oxygen of phosphate, called 1BPh
    % 11  G N2-2H2 interacts with oxygen of phosphate, called 3BPh
    % 12  G N2-2H2 and N1-H1 interacts with 2 oxygens of phosphate, called 4BPh
    % 19  G N2-2H2 and N1-H1 interact with just one oxygen, called 4BPh
    % 13  G N1-H1  interacts with oxygen of phosphate, called 5BPh
    % 14  G C8-H8  interacts with oxygen of phosphate, called 0BPh
    %
    % 15  U N3-H3  interacts with oxygen of phosphate, called 5BPh
    % 16  U C5-H5  interacts with oxygen of phosphate, called 9BPh
    % 17  U C6-H6  interacts with oxygen of phosphate, called 0BPh
"""
    CONDITIONS = {
        'A': {
            '2': [[('i_C2_H2',1)]],
            '6': [[('i_N6_1H6',1)]],
            '7': [[('i_N6_2H6',1)]],
            '789': [[('i_N6_2H6',1)]],
            '0': [[('i_C8_H8',1)]],
        },
        'C': {
            '6': [[('i_N4_2H4',1)]],
            '789': [
                    [('i_N4_1H4',1)],
                    [('i_C5_H5',1)],
                 ],
            '8': [[('i_N4_1H4',1),('i_C5_H5',1),('oxygens_count',2,100)]],
            '7': [
                    [('i_N4_1H4',1),('i_C5_H5',0)],
                    [('i_N4_1H4',1),('i_C5_H5',1),('oxygens_count',1)]
                 ],
            '9': [[('i_N4_1H4',0),('i_C5_H5',1),('i_C6_H6',0,1)]],
            '0': [[               ('i_C5_H5',0),('i_C6_H6',1)]],
        },
        'G': {
            '1': [[('i_N2_1H2',1)]],
            '345': [
                    [('i_N2_2H2',1)],
                    [('i_N1_H1',1)],
                 ],
            '3': [[('i_N2_2H2',1),('i_N1_H1',0)]],
            '4': [[('i_N2_2H2',1),('i_N1_H1',1),('oxygens_count',1,100)]],
            '5': [[('i_N2_2H2',0),('i_N1_H1',1)]],
            '0': [[('i_C8_H8',1)]],
        },
        'U': {
            '5': [[('i_N3_H3',1)]],
            '345': [[('i_N3_H3',1)]],
            '9': [[('i_C5_H5',1),('i_C6_H6',0,1)]],
            '789': [[('i_C5_H5',1),('i_C6_H6',0,1)]],
            '0': [[('i_C5_H5',0),('i_C6_H6',1)]],
        },
    }

    if options.sub_category=='base-phosphate':
        ox_set = PH_OXYGENS
    elif options.sub_category=='base-ribose':
        ox_set = BR_OXYGENS
    else:
        raise Exception("Unknown sub-category=%s"%options.sub_category)

    prefix=""
    if options.sub_category=='base-ribose':
        prefix="br_"
    elif options.sub_category=='base-phosphate':
        prefix="ph_"

    
    # fuzzy detection
    (ref,fuzzy,_other,unclassified) = load_doublets(options.sub_category,options.desc,options.n_type,options)
    print "loaded %s reference doublets, %s unclassified" % (len(ref),len(unclassified))
    if len(unclassified)>10000:
        print "truncating unclassified set to 10000"
        unclassified = list(unclassified)
        random.shuffle(unclassified)
        unclassified = unclassified[0:10000]
    ids = set(list(ref)+list(unclassified))
    dd = DoubletsDict(options.data_dir,reduced_atoms=R_ATOMS)
    dd.load_pdb_files(ids,verbose=True)
    all_params = compute_doublet_params(ref, options.sub_category, options)

    if len(ref)>0:
        cl_doublets = [x for x in ref]
        random.shuffle(cl_doublets)
        cl_doublets = cl_doublets[0:8]
        save_doublets_as_pdb(fn("doublets-pdb",sc=options.sub_category,desc=options.desc,n_type=options.n_type,opt="cl"), cl_doublets, "full")
        save_json(fn("doublets-json",sc=options.sub_category,desc=options.desc,n_type=options.n_type,opt="cl"), cl_doublets)

    points = []
    base = options.n_type[0].upper()
    num = re.sub('[^0-9]','',options.desc)
    for d_id in ref:
        p = dd.get_normalized(d_id)
        params = all_params.get(d_id)
        if params is None:
            print "WARNING! skipping %s, missing params" % d_id
            continue
        for C in CONDITIONS[base][num]:
            for c in C:
                if len(c)==2 and c[1]==1:
                    for k,v in params.items():
                        if re.match('^i'+c[0]+'_(.*)',k) and v==1:
                            atom_name = k.split("_")[3]
                            points.append(p[1][atom_name])
    print "ref_points: %s" % len(points)


    # exact detection
    d_limits_exact = None
    if len(points)>0:
        d_limits_exact = point_classifier(
                array(points,'f'),
                array([(0,0,0)],'f'),
                array([(0,0,0)],'f'),
                outliers_count=0,
                cl_count=50,
                cl_groups=1,
                dist_factor=2.0
        )
    
    ret_code = []
    if d_limits_exact is not None:
        for o_atom in ox_set:
            for limit_num,(d,p) in enumerate(d_limits_exact,start=1):
                ret_code += [["COMPUTE_DISTANCE_TO_POINTS","$norm2_%s"%o_atom,p]]
                ret_code += [["RETURN_IF","$distance_to_points",d,d+0.00001,"%s (%s-%d)"%(options.desc,o_atom,limit_num)]]
    else:
        ret_code += [["RETURN",options.desc]]

    base = options.n_type[0].upper()
    num = re.sub('[^0-9]','',options.desc)
    test_code = []
    if CONDITIONS.has_key(base) and CONDITIONS[base].has_key(num):
        for C in CONDITIONS[base][num]:
            cond = []
            cond.append(["EQ","$n_type",options.n_type])
            if options.sub_category in ['base-phosphate']:
                cond.append(["IN_RANGE","$d_param_%s%s"%(prefix,"ph_h"),0,4.5])
            for c in C:
                if len(c)==2:
                    cond.append(["EQ","$d_param_%s%s"%(prefix,c[0]),c[1]])
                elif len(c)==3:
                    cond.append(["IN_RANGE","$d_param_%s%s"%(prefix,c[0]),c[1],[2]])
            test_code.append(["IF",cond,ret_code])


    # fuzzy detection
    if len(points)>0:
        uncl_points = []
        for d_id in unclassified:
            p = dd.get_normalized(d_id)
            for o_atom in ox_set:
                if p[1].has_key(o_atom):
                    uncl_points.append(p[1][o_atom])
        print "uncl_points: %s" % len(uncl_points)
                    
        d_limits = point_classifier(
                    array(points,'f'),
                    array([(0,0,0)],'f'),
                    array(uncl_points,'f'),
                    outliers_count=500,
                    cl_count=800,
                    cl_groups=5,
                    dist_factor=1.2
        )
        uncl_points = filter_points_by_distance(uncl_points, d_limits) 
        print "after distance filtering, uncl_points: %s" % len(uncl_points)
    
        fuzzy_cl = []    
        for o_atom in ox_set:
            for num,(d,p) in enumerate(d_limits,start=1):
                fuzzy_cl += [["COMPUTE_DISTANCE_TO_POINTS","$norm2_%s"%o_atom,p]]
                fuzzy_cl += [["RETURN_IF","$distance_to_points",-0.001,d,"%s (%s-%d)"%(options.desc,o_atom,num)]]
        test_code += [["IF",[["EQ","$n_type",options.n_type]],fuzzy_cl]]
                            
        save_points_as_pdb(fn("doublets-pdb",sc=options.sub_category,desc=options.desc,n_type=options.n_type,opt="ref-oxygens"), points, options.n_type)
        save_points_as_pdb(fn("doublets-pdb",sc=options.sub_category,desc=options.desc,n_type=options.n_type,opt="uncl-oxygens"), uncl_points, options.n_type)

    save_doublets_as_pdb(fn("doublets-pdb",sc=options.sub_category,desc=options.desc,n_type=options.n_type,opt="ref"), ref, options.sub_category)
    
    f = fn("classifier",sc=options.sub_category,n_type=options.n_type,desc=options.desc)
    print "saving output to %s" % f
    save_json(f, test_code, indent=True)
    print "OK"

def compute_other_spec_classifier(options):
    (ref,fuzzy,other,unclassified) = load_doublets(options.sub_category,options.desc,options.n_type,options)
    if options.desc in ['phosphate-stacking1','phosphate-stacking2']:
        unclassified = set(list(unclassified)[0:10])
    print "loading doublets, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified))

    if options.desc in ['long-stacking-c','diagonal-c','long-diagonal-c','diagonal-nc-ww']:
        def only_consecutive(elems,exp_value=1):
            res = []
            for x in elems:
                _pdbid,x1,x2 = x.split(":")
                if x1[0]==x2[0] and abs(int(x1[1:])-int(x2[1:]))==1:
                    if exp_value==1:
                        res.append(x)
                else:
                    if exp_value==0:
                        res.append(x)
            return res
        exp_value = 1
        if options.desc in ['diagonal-nc-ww']:
            exp_value = 0
        ref = set(only_consecutive(ref, exp_value))
        fuzzy = set(only_consecutive(fuzzy, exp_value))
        other = set(only_consecutive(other, exp_value))
        unclassified = set(only_consecutive(unclassified, exp_value))
        print "after consecutive filter, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified))
    unclassified = set(list(unclassified)[0:1000])

    all_ids = ref.union(fuzzy).union(unclassified)
    param_type = 'bp'
    param_set = ['nn_ang_norm','dist','n12cc_ang']
    if options.desc in ['base-ribose-stacking']:
        param_type = 'base-ribose'
        param_set = ['br_dist','n1br_ang_norm']
    elif options.desc in ['phosphate-stacking1','phosphate-stacking2']:
        param_type = 'base-phosphate'
        param_set = ['n1ph_ang']
    all_params = compute_doublet_params(all_ids, param_type, options)
    print all_params
        
    for k in param_set:
        p = [all_params[d_id][k] for d_id in ref]
        if len(p)==0:
            print "empty set of key: %s" % k
            continue
        if len(p)>100:
            n = len(p)
            p.sort()
            p = p[max(n//100,1):n-max(n//100,1)]
        print k, min(p), max(p)

    limits = None
    limit_consecutive = None
    param_prefix = ""
    if options.desc=='long-stacking-c':
        limit_consecutive = True
        limits = {'nn_ang_norm': {'min':0.0,'max':40.0}, 
                  'dist': {'min':6.4,'max':9.0}, 
                  'n12cc_ang': {'min': 0.0,'max':45.0}, 
                  }
    elif options.desc=='diagonal-c':
        limit_consecutive = True
        limits = {'nn_ang_norm': {'min':0.0,'max':60.0}, 
                  'dist': {'min':4.0,'max':9.0}, 
                  'n12cc_ang': {'min': 30.0,'max':70.0}, 
                  }
    elif options.desc=='long-diagonal-c':
        limit_consecutive = True
        limits = {'nn_ang_norm': {'min':0.0,'max':41.0}, 
                  'dist': {'min':8.0,'max':12.0}, 
                  'dist_z': {'min':4.7,'max':8.0}, 
                  'n12cc_ang': {'min': 40.0,'max':60.0}, 
                  }
    elif options.desc=='diagonal-nc-ww':
        limit_consecutive = False
        limits = {'nn_ang_norm': {'min':0.0,'max':40.0}, 
                  'dist': {'min':5.5,'max':10.5}, 
                  'dist_z': {'min':1.8,'max':5.5}, 
                  'n12cc_ang': {'min': 50.0,'max':80.0}, 
                  }
    elif options.desc=='base-ribose-stacking':
        param_prefix = "br_"
        limits = {'br_dist': {'min':3.0,'max':6.0}, 
                  'n1br_ang_norm': {'min':0.0,'max':40.0}, 
                  }
    elif options.desc=='phosphate-stacking1':
        param_prefix = "ph_"
        limits = {'n1ph_ang': {'min':0.0,'max':30.0},'ph_dist': {'min':0.0,'max':5.0}}
    elif options.desc=='phosphate-stacking2':
        param_prefix = "ph_"
        limits = {'n1ph_ang': {'min':150.0,'max':180.0},'ph_dist': {'min':0.0,'max':5.0}}

    
    if limits is not None:
        org_ref = ref.copy()
        ref = filter_by_params(ref, limits, all_params) 
        other = filter_by_params(other, limits, all_params) 
        fuzzy = filter_by_params(fuzzy, limits, all_params) 
        unclassified = filter_by_params(unclassified, limits, all_params) 
    
        print "after filtering by params, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified))
        if len(org_ref)!=len(ref):
            removed = set(org_ref).difference(set(ref))
            print "removed doublets(%d): %s" % (len(removed),list(removed)[0:10])
            if False:
                expert_fn = "expert/ref/other_%s_%s.json" % (options.n_type,options.desc)
                data = load_json(expert_fn)
                data = [x for x in data if x not in removed]
                print "new data=%s" % len(data)
                save_json(expert_fn,data)

    d_limits = None
    if len(ref)>0 and options.desc in ['long-stacking-c','diagonal-c','long-diagonal-c','diagonal-nc-ww']:
    
        cl_count=500
        dist_factor=1.2
        
        if options.desc in ['diagonal-c'] and options.n_type in ['GU']:
            cl_count=800
            dist_factor=1.6
        elif options.desc in ['diagonal-c'] and options.n_type in ['UC','CU']:
            cl_count=800
            dist_factor=1.8
        
        d_limits = point_classifier(
            array([all_params[d_id]['center'] for d_id in ref],'f'),
            array([(0,0,0)],'f'),
            array([all_params[d_id]['center'] for d_id in unclassified],'f'),
            outliers_count=1000,
            cl_count=cl_count,
            cl_groups=5,
            dist_factor=dist_factor
        )

    if d_limits is not None:
        ref = filter_by_distance(ref, d_limits, all_params) 
        other = filter_by_distance(other, d_limits, all_params) 
        fuzzy = filter_by_distance(fuzzy, d_limits, all_params) 
        unclassified = filter_by_distance(unclassified, d_limits, all_params) 
        print "after filtering by distance, ref=%d fuzzy=%d other=%d unclassified=%d" % (len(ref), len(fuzzy), len(other), len(unclassified))

    cl = [["RETURN",options.desc]]

    if d_limits is not None:
        cl = []
        for num,(d,p) in enumerate(d_limits,start=1):
            cl += [["COMPUTE_DISTANCE_TO_POINTS","$norm2_center",p]]
            cl += [["RETURN_IF","$distance_to_points",d,4*d,"%s (%d)"%(options.desc,num)]]

    if limits is not None:
        cond = []
        if limit_consecutive is not None:
            if limit_consecutive:
                cond += [["EQ","$is_consecutive",1]]
            else:
                cond += [["EQ","$is_consecutive",0]]
        cond += [["IN_RANGE","$d_param_%s%s"%(param_prefix,k),v['min'],v['max']] for k,v in limits.items()]
        cl = [["IF",cond,cl]]

    test_code = [["IF",[["EQ","$n_type",options.n_type]],cl]]
                            
    # save_doublets_as_pdb(fn("doublets-pdb",sc=options.sub_category,desc=options.desc,n_type=options.n_type,opt="ref"), ref, options.sub_category)
    
    f = fn("classifier",sc=options.sub_category,n_type=options.n_type,desc=options.desc)
    print "saving output to %s" % f
    save_json(f, test_code, indent=True)
    print "OK"


def test_classifier(options):
    print "test!"
    cl_fn = fn("classifier",sc=options.sub_category,n_type=options.n_type,desc=options.desc)
    cl_lib = __import__("clarna")
    cl = cl_lib.ContactProgram()
    cl.load_from_json(cl_fn)

    (ref,fuzzy,other,unclassified) = load_doublets(options.sub_category,options.desc,options.n_type,options)
    
    all_ids = ref.union(fuzzy).union(unclassified)
    all_params = compute_doublet_params(all_ids, options.sub_category, options)
    
    class Opt: pass
    cl_opt = Opt()
    
    results = dict([(k1+"_"+k2,0) for k1 in ['ref','fuzzy','other','unclassified'] for k2 in ['detected','undetected','fuzzy']])
    
    i = 0
    def show_results(r,stderr=False):
        res = "RESULTS:"
        for k in ['ref','fuzzy','other','unclassified']:
            res += " "+k+"=%d,%d,%d" % (r[k+"_detected"], r[k+"_fuzzy"], r[k+"_undetected"])
        if stderr:
            print >>sys.stderr, res
        else:
            print res
    
    for name,doublets in ('ref',ref),('fuzzy',fuzzy),('other',other),('unclassified',unclassified):
        for d_id in doublets:
            if not all_params.has_key(d_id):
                continue
            d_params = all_params[d_id]
            normalized_r1 = {}
            normalized_r2 = {"center": d_params.get('norm2_center')}
            state = cl_lib.EmptyContactProgramState(n_type=options.n_type, doublet_params=d_params, normalized_r2=normalized_r2)
            (res,_state) = cl.value(state)
            res = cl_lib.single_result(res)
            if res is not None and res.desc[0]!='?':
                results[name+"_detected"] += 1
            elif res is not None and res.desc[0]=='?':
                results[name+"_fuzzy"] += 1
            else:
                results[name+"_undetected"] += 1
            i += 1
            if i%500==0:
                show_results(results,stderr=True)
    show_results(results)

def main():
    random.seed(12345)

    (parser,options,_args) = parse_args()
    
    global OUTPUT_DIR
    if options.output_dir:
        OUTPUT_DIR = options.output_dir
        
    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
    
    # test1('CA',options)
    # test2('WW','CA',options)
    # test3('CA',options)
    
    #options.n_type = 'CA'
    #options.sub_category = 'bp'

    #options.desc = "WW"
    
    if options.combine:
        return combine_classifier(options)
        
    if options.test_classifier:
        return test_classifier(options)

    
    assert options.n_type is not None
    assert options.sub_category is not None
    
    if options.desc is None:
        if options.sub_category=='bp':
            compute_bp_classifier(options)
        else:
            raise Exception("Unsupported subcategory: %s" % options.sub_category)
    else:
        if options.sub_category=='bp':
            assert re.match('^[SWH]{2}(_cis|_tran)?$',options.desc)
            compute_bp_spec_classifier(options)
        elif options.sub_category=='stacking':
            compute_stacking_spec_classifier(options)
        elif options.sub_category in ['base-phosphate','base-ribose']:
            bph_and_br_fr3d_classifier(options)
        elif options.sub_category in ['other','other2','other3']:
            compute_other_spec_classifier(options)
        else:
            raise Exception("Unsupported: %s/%s" % (options.sub_category,options.desc))
        
    
if __name__=="__main__":
    main()