#!/usr/bin/env python
import sys
import os
import re
import math
import gzip
import itertools
from optparse import OptionParser
from itertools import combinations
import shutil
import networkx as nx
from networkx.readwrite import json_graph
import numpy as np

from utils import load_json,save_json, GroupsTool, GraphTool, DoubletDescTool, DoubletsDict, confusion_matrix_params, status_msg
from distances import doublet_params_dict, expected_strand_orient

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""generate the json dictionaries 
with comparison of classifier graph to the reference groups""")
    parser.add_option("--groups", dest="groups", metavar="JSON_FILE")
    parser.add_option("--classifier-graph", dest="classifier_graph", metavar="FILE")
    parser.add_option("--close-doublets-graph", dest="close_doublets_graph", metavar="FILE")
    parser.add_option("--rnaview-graph", dest="rnaview_graph", metavar="FILE")
    parser.add_option("--mc-annotate-graph", dest="mc_annotate_graph", metavar="FILE")
    parser.add_option("--moderna-graph", dest="moderna_graph", metavar="FILE")
    parser.add_option("--fr3d-graph", dest="fr3d_graph", metavar="FILE")
    parser.add_option("--pdb-id", dest="pdb_id", metavar="ID")
    parser.add_option("--eval-mode", dest="eval_mode", action='store_true', default=False)
    parser.add_option("--new-eval-mode", dest="new_eval_mode", action='store_true', default=False)
    parser.add_option("--correlation", dest="correlation", action='store_true', default=False)
    parser.add_option("--verbosity", dest="verbosity", default=0, metavar="N")
    parser.add_option("--only-doublets-from", dest="only_doublets_from", metavar="FILE")
    parser.add_option("--data-dir", dest="data_dir",
                  help="directory with data", metavar="DIR", default="gc-data")
    parser.add_option("-o", "--output-json", dest="output_json",
                  help="save result to output JSON", metavar="FILE")

    (options, args)  = parser.parse_args()
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

def load_doublets_from_file(options):
    pdb_id = options.pdb_id.upper()
    fn = options.only_doublets_from
    data = load_json(fn)
    res = set()
    for elems in data.values():
        for x in elems:
            xx = x.split(":")
            if xx[0]==pdb_id:
                res.add(xx[1]+":"+xx[2])
                res.add(xx[2]+":"+xx[1])
    return res

def graph_to_doublets(g,prg=None):
    doublets = {}
    if g is None:
        return doublets
        
    # temporary fix for multiple labels for single contact (for example in FR3D)
    def _edge_key(source,target,data):
        desc = data.get('desc','')
        desc_key = "0"
        if re.match("^[sSWH]{2}_(cis|tran)$", desc):
            desc_key = "9"
        elif re.match("^[<>]{2}$", desc):
            desc_key = "5"
        elif re.match("^[SWH]{1,2}_[0-9]+", desc):
            desc_key = "2"
        return "%s:%s:%s:%s" % (source,target,desc_key,desc)
    
    for (source,target,data) in sorted(g.edges(data=True), key=lambda x: _edge_key(x[0],x[1],x[2])):
        if data.get('reverse',False)==False:
            id="%s:%s"%(source,target)
            if prg is not None:
                if data['type']=='contact' and data['prg']!=prg:
                    continue
            doublets[id]={'desc':data.get('desc'),'desc_full':data.get('full_desc'),'type':data['type'],'n_type':data.get('n_type')}
    return doublets

def get_expected_result(rv_row,mc_row,mo_row,fr_row):
    """TODO: use desc tool"""
    def _get(x):
        if x is not None:
            if x['desc'] is not None:
                return x['desc']
            else:
                return ''
        else:
            return ''
    rv_res = _get(rv_row)
    mc_res = _get(mc_row)
    mo_res = _get(mo_row)
    fr_res = _get(fr_row)
    if rv_res==mc_res and re.match("^[sHSW]{2}_(cis|tran)$",rv_res):
        return rv_res
    elif mc_res==mo_res and re.match("^[<>]{2}$",mc_res):
        return mc_res
    elif re.match("^[A-Z]+_[0-9]",fr_res):
        return fr_res
    if rv_res!='' or mc_res!='' or mo_res!='' or fr_res!='':
        a = []
        for res in [rv_res,mc_res,mo_res,fr_res]:
            if res!='' and res!='UNK_SHORT_RES':
                a.append(res)
        return '?('+("|".join(a))+')'
    return ''

def add_result(results,sc,_id,desc,n_type):
    desc = desc.replace('?','')
    key = 'evaluation/'+sc+'/'+desc+'/'+n_type
    if not results.has_key(key):
        results[key]=[]
    # print "adding %s to %s" % (_id,key)
    results[key].append(_id)

def combine_graphs_old(options):
    cl_graph = load_graph(options.classifier_graph)
    ###
    rv_graph = load_graph(options.rnaview_graph)
    mc_graph = load_graph(options.mc_annotate_graph)
    mo_graph = load_graph(options.moderna_graph)
    fr_graph = load_graph(options.fr3d_graph)
    ###
    cl_d = graph_to_doublets(cl_graph,'MY')
    rv_d = graph_to_doublets(rv_graph)
    mc_d = graph_to_doublets(mc_graph)
    mo_d = graph_to_doublets(mo_graph)
    fr_d = graph_to_doublets(fr_graph)
    ###
    results = {}
    for _id,d in cl_d.items():
        desc = d.get('desc')
        if desc is None:
            desc=''
        n_type = d.get('n_type')
        exp_res = get_expected_result(rv_d.get(_id),mc_d.get(_id),mo_d.get(_id),fr_d.get(_id))
        full_id = _id
        if options.pdb_id:
            full_id = options.pdb_id.upper() + ":" + _id
        else:
            full_id = _id
        if exp_res!='' and exp_res[0]!='?':
            # print "reference doublet %s: %s" % (_id,exp_res)
            add_result(results,'ref-all',full_id,exp_res,n_type)
            if exp_res==desc:
                add_result(results,'ref-ok',full_id,exp_res,n_type)
            elif ('?'+exp_res)==d['desc']:
                add_result(results,'ref-ok-fuzzy',full_id,exp_res,n_type)
            elif desc=='':
                add_result(results,'ref-undetected',full_id,exp_res,n_type)
            else:
                add_result(results,'ref-diff',full_id,exp_res,n_type)
        elif exp_res!='' and exp_res[0]=='?' and desc!='' and desc[0]=='?':
            exp_results = exp_res[2:len(exp_res)-1].split("|")
            # print desc[1:], exp_results
            if desc[1:] in exp_results:
                add_result(results,'fuzzy-ok',full_id,desc[1:],n_type)
        elif exp_res=='' and desc!='' and desc[0]!='?':
            add_result(results,'prev-undetected',full_id,desc,n_type)
    save_json(options.output_json, results)

def combine_graphs_using_groups(options):
    filter_doublets = None
    if options.only_doublets_from:
        filter_doublets = load_doublets_from_file(options)

    groups = GroupsTool(options.groups,filter_doublets=filter_doublets)
    cl_graph = GraphTool(options.classifier_graph,filter_doublets=filter_doublets)
    ###
    results = {}
    
    all_ids = list(groups.get_all_ids())
    recognized_by = dict((prg,groups.get_recognized_by(prg)) for prg in ['RV','MC','MO','FR'])
    for full_id in sorted(all_ids):
        n_type = groups.get_n_type(full_id)
        short_id = re.sub("^[^:]+:","", full_id)
        rev_full_id = DoubletDescTool.reverse_d_id(full_id)
        rev_short_id = DoubletDescTool.reverse_d_id(short_id)
        cl_desc = cl_graph.get_contact_by_id(short_id)
        rev_cl_desc = cl_graph.get_contact_by_id(rev_short_id)
        recognized_by_cl = False
        if cl_desc!='' and rev_cl_desc!='':
            recognized_by_cl = True
        if not recognized_by_cl:
            add_result(results,'unclassified',full_id,'by-classifier',n_type)

        recognized_by_any = recognized_by_cl
        for prg in ['RV','MC','MO','FR']:
            if full_id in recognized_by[prg] or rev_full_id in recognized_by[prg]:
                recognized_by_any = True
            else:
                # add_result(results,'unclassified',full_id,'by-%s'%prg.lower(),n_type)
                pass
        if not recognized_by_any:
            add_result(results,'unclassified',full_id,'by-any',n_type)
    
    for sub_category in ('bp','stacking','base-phosphate','base-ribose','other','other2','other3'):
        ref_ids = set(groups.get_ref_ids(sub_category))
        fuzzy_ids = set(groups.get_fuzzy_ids(sub_category))

        for full_id in ref_ids:
            n_type = groups.get_n_type(full_id)
            exp_desc = groups.get_ref_desc(full_id, sub_category)
            assert exp_desc != ''
            assert n_type != '??'

            if sub_category in ['bp','stacking']:
                rev_n_type = DoubletDescTool.reverse_n_type(n_type)
                rev_exp_desc = DoubletDescTool.reverse_desc(exp_desc)
                if n_type==rev_n_type and exp_desc==rev_exp_desc and GroupsTool.is_reversed_id(full_id):
                    continue


            add_result(results,'ref-all',full_id,exp_desc,n_type)

            short_id = re.sub("^[^:]+:","", full_id)
            cl_desc = cl_graph.get_contact_by_id(short_id, sub_category)

            if cl_desc == exp_desc:
                add_result(results,'ref-ok',full_id,exp_desc,n_type)
            elif cl_desc!='' and cl_desc[0]=='?' and cl_desc[1:]==exp_desc:
                add_result(results,'ref-ok-fuzzy',full_id,exp_desc,n_type)
            elif cl_desc == '' or cl_desc[0]=='?':
                add_result(results,'ref-undetected',full_id,exp_desc,n_type)
            else:
                add_result(results,'ref-diff',full_id,exp_desc,n_type)

        for full_id in fuzzy_ids:
            if GroupsTool.is_reversed_id(full_id):
                continue
            n_type = groups.get_n_type(full_id)
            exp_desc = groups.get_fuzzy_desc(full_id, sub_category)
            assert len(exp_desc)!=0
            if n_type=='??':
                print "ERROR!", full_id, exp_desc, n_type
            assert n_type != '??'

            short_id = re.sub("^[^:]+:","", full_id)
            cl_desc = cl_graph.get_contact_by_id(short_id, sub_category)

            if cl_desc in exp_desc:
                add_result(results,'fuzzy-ok',full_id,cl_desc,n_type)

        for short_id in cl_graph.get_ids(sub_category):
            full_id = options.pdb_id.upper()+":"+short_id

            cl_desc_dict = cl_graph.get_contact_by_id(short_id,sub_category,data=True)
            cl_desc = cl_desc_dict['desc']
            assert cl_desc != ''

            if cl_desc[0]=='?':
                # add_result(results,'detected-by-cl-fuzzy',full_id,cl_desc,cl_desc_dict['n_type'])
                pass
            else:
                # add_result(results,'detected-by-cl',full_id,cl_desc,cl_desc_dict['n_type'])
                pass
            

            if full_id in ref_ids or full_id in fuzzy_ids:
                continue
            
            if cl_desc[0]=='?':
                # add_result(results,'prev-undetected-fuzzy',full_id,cl_desc[1:],cl_desc_dict['n_type'])
                pass
            else:
                add_result(results,'prev-undetected',full_id,cl_desc,cl_desc_dict['n_type'])

    save_json(options.output_json, results)

# new code!!!!
def _combine_graphs_eval_stat2(results, graphs, sub_categories, options):
    pdb_id = options.pdb_id.upper()
    all_programs = sorted(graphs.keys())
    verbosity = int(options.verbosity)
    for sc in sub_categories:
        ids = {}
        id_data = {}
        for prg in all_programs:
            if sc in ['bp','stacking']:
                ids[prg] = set(graphs[prg].get_unique_ids(sc))
            else:
                ids[prg] = set(graphs[prg].get_ids(sc))
            
            id_data[prg] = {}
            for id in ids[prg]:
                id_data[prg][id] = graphs[prg].get_contact_by_id(id,sc,data=True)
        descriptions = set([d['desc'] for prg in all_programs for d in id_data[prg].values()])
        for n in xrange(1,len(all_programs)+1):
            for sel_programs in combinations(all_programs, n):
                key1 = "".join(sel_programs)
                prg0 = sel_programs[0]
                for id in ids[prg0]:
                    ok = True
                    desc = id_data[prg0][id]['desc']
                    n_type = id_data[prg0][id]['n_type']
                    if not re.match("^[ACGU]{2}$",n_type):
                        continue
                    for prg in all_programs:
                        if prg in sel_programs:
                            if (not id in ids[prg]) or (desc!=id_data[prg][id]['desc']):
                                ok = False
                                break
                        else: # prg not in sel_programs
                            if (id in ids[prg]) and (desc==id_data[prg][id]['desc']):
                                ok = False
                                break
                    if ok:
                        key = "stat2/"+key1+"/"+sc+"/"+desc+"/"+n_type
                        if not results.has_key(key):
                            results[key]=[]
                        results[key].append(pdb_id+":"+id)
    if verbosity>0:
        for k,v in results.items():
            if k[0:5]=="stat2":
                print "stat2",k,len(v)

def combine_graphs_eval(options):
    filter_doublets = None
    if options.only_doublets_from:
        filter_doublets = load_doublets_from_file(options)

    g_cl = GraphTool(options.classifier_graph, accept_fuzzy=False, filter_doublets=filter_doublets)
    g_mc = GraphTool(options.mc_annotate_graph, filter_doublets=filter_doublets)
    g_fr = GraphTool(options.fr3d_graph, filter_doublets=filter_doublets)
    g_mo = GraphTool(options.moderna_graph, filter_doublets=filter_doublets)

    graphs_rv = GraphTool()
    graphs_rv.load_graph(options.rnaview_graph)
    for v1,v2,data in graphs_rv.g.edges(data=True):
        if data.has_key('full_desc') and "stacking" in data['full_desc']:
            data['desc']='<<'
    graphs_rv.contacts = GraphTool._graph_to_contacts(graphs_rv.g)

    graphs = {
        "CL": g_cl,
        "MC": g_mc,
        "FR": g_fr,
        "RV": graphs_rv,
    }

    
    graphs_st = {
        "CL": g_cl,
        "MC": g_mc,
        "FR": g_fr,
        "MO": g_mo,
    }
    
    sub_categories = ('bp','stacking','base-phosphate','base-ribose')
    results = {}
    _combine_graphs_eval_stat2(results, graphs, ('bp','base-phosphate','base-ribose'), options)
    _combine_graphs_eval_stat2(results, graphs_st, ('stacking',), options)
    save_json(options.output_json, results)

def combine_graphs_correlation(options):
    filter_doublets = None
    if options.only_doublets_from:
        filter_doublets = load_doublets_from_file(options)

    verbosity = int(options.verbosity)

    CATEGORIES = ["bp-classic","bp-non-classic","stacking","base-phosphate","base-ribose"]
    PRG_BY_CAT = {
        "bp-classic": ["CL","RV","MC","FR"],
        "bp-non-classic": ["CL","RV","MC","FR"],
        "stacking": ["CL","MC","MO","FR"],
        "base-phosphate": ["CL","FR"],
        "base-ribose": ["CL","FR"],
    }

    graphs = {}
    for prg_code,fn in [
            ('CL',options.classifier_graph),
            ('FR',options.fr3d_graph),
            ('MC',options.mc_annotate_graph),
            ('MO',options.moderna_graph),
            ('RV',options.rnaview_graph),
        ]:
        graphs[prg_code] = GraphTool(fn,accept_fuzzy=False,filter_doublets=filter_doublets)
    
    res = {}
    for sc in CATEGORIES:
        P_LIST = PRG_BY_CAT[sc]

        for (p1,p2) in itertools.product(sorted(P_LIST), repeat=2):
            if sc in ['bp-classic','bp-non-classic','stacking']:
                p1_ids = set(graphs[p1].get_unique_ids(sc))
                p2_ids = set(graphs[p2].get_unique_ids(sc))
            else:
                p1_ids = set(graphs[p1].get_ids(sc))
                p2_ids = set(graphs[p2].get_ids(sc))
            common_ids = []
            for d_id in p1_ids.intersection(p2_ids):
                if graphs[p1].get_contact_by_id(d_id,cat=sc)==graphs[p2].get_contact_by_id(d_id,cat=sc):
                    common_ids.append(d_id)
            all_ids = p1_ids.union(p2_ids)
            key1 = "correlation/"+sc+"/"+p1+"_vs_"+p2+"/common"
            key2 = "correlation/"+sc+"/"+p1+"_vs_"+p2+"/all"

            if not res.has_key(key1):
                res[key1] = 0
            res[key1] += len(common_ids)
            if not res.has_key(key2):
                res[key2] = 0
            res[key2] += len(all_ids)
            if verbosity>0:
                status_msg("correlation/%s/%s_vs_%s = %.4f (common=%d,all=%d)" % (sc,p1,p2,float(len(common_ids))/max(1,len(all_ids)), len(common_ids), len(all_ids)))
            
    save_json(options.output_json, res)

class GraphsCollection:

    def __init__(self,pdb_id,cl_fn,fr_fn,mc_fn,mo_fn,rv_fn,dd=None,filter_doublets=None):
        self.pdb_id = pdb_id.lower()
        self.graphs = {}
        self.graphs["CL"] = GraphTool(cl_fn,accept_fuzzy=True,filter_doublets=filter_doublets)
        for prg_code,prg_fn in [('FR',fr_fn),('MC',mc_fn),('MO',mo_fn),('RV',rv_fn)]:
            self.graphs[prg_code] = GraphTool(prg_fn,accept_fuzzy=False,filter_doublets=filter_doublets)
        self.dd = dd
        self.PRG_BY_CAT = {
            "bp": ["CL","RV","MC","FR"],
            "bp-classic": ["CL","RV","MC","FR"],
            "bp-non-classic": ["CL","RV","MC","FR"],
            "stacking": ["CL","MC","MO","FR"],
            "base-phosphate": ["CL","FR"],
            "base-ribose": ["CL","FR"],
        }

    def _strand_orientation_name(self,x):
        if x==1:
            return "parallel"
        elif x==-1:
            return "anti-parallel"
        else:
            raise Exception("bad value: %s" % x)


    def get_contact_by_id(self,prg,short_id,cat,data=False):
        return self.graphs[prg].get_contact_by_id(short_id,cat,data)
    
    def classify_prg_result_details_fp(self,short_id,sc,prg_res,exp_res,other_results):
        other_prgs = [x for x in self.PRG_BY_CAT[sc] if x!='CL']

        d_id = self.pdb_id.upper()+":"+short_id
        if sc=='base-ribose':
            p_br = doublet_params_dict(self.dd.get(d_id), self.dd.get_n_type(d_id), 'base-ribose')
        if sc in ['bp-classic','bp-non-classic']:
            p_bp = doublet_params_dict(self.dd.get(d_id), self.dd.get_n_type(d_id), 'bp')
        
        sc2 = sc.replace("bp-classic","bp").replace("bp-non-classic","bp")
        
        exp_res2 = exp_res
        if exp_res2 == "":
            oo = [x for x in other_results if x!=""]
            if len(oo)==1:
                exp_res2 = oo[0]

        if sc=='base-ribose' and p_br['ph_h']>=4.5:
            return ("new-base-ribose","ph_h=%.4f"%p_br['ph_h'])
        elif sc2=='bp' and prg_res!="" and exp_res2!="" and prg_res[0:2]==exp_res2[0:2] and expected_strand_orient(exp_res2)!=p_bp['strand_orient']:
            extra = "strand orientation: %s, expected orientation for %s: %s"% (
                        self._strand_orientation_name(p_bp['strand_orient']),
                        exp_res2,
                        self._strand_orientation_name(expected_strand_orient(exp_res2))
                    )
            return ("cis-vs-trans",extra)
        elif sc2=='bp' and prg_res=='SH_cis' and exp_res2=='WH_cis':
            extra = ""
            return ("wh_cis-vs-sh_cis",extra)
        else:
            if prg_res in other_results:
                extra = "consistent with: "+[p for p,o in zip(other_prgs,other_results) if prg_res==o][0]
                return ("consistent-with-single-cl",extra)
            elif all([x=="" for x in other_results]):
                return ("not-recognized-by-others-cl","")
            else:
                return ("others","")
    
    def classify_prg_result(self,prg,short_id,sc,details=False,min_contact_weight=None):
        assert self.graphs.has_key(prg)
        assert self.PRG_BY_CAT.has_key(sc)

        d_id = self.pdb_id.upper()+":"+short_id
        prg_res = self.get_contact_by_id(prg,short_id,cat=sc)
        prg_res_data = self.get_contact_by_id(prg,short_id,cat=sc,data=True)
        if min_contact_weight is None:
            if "?" in prg_res:
                prg_res = ""
        else:
            prg_res = prg_res.replace("?","")
            w = prg_res_data.get('weight')
            if w is None:
                if "?" in prg_res:
                    w = 0.4999
                else:
                    w = 1.0
            if w<min_contact_weight:
                prg_res = ""
        
        other_prgs = [x for x in self.PRG_BY_CAT[sc] if x!=prg and x!='CL']
        other_results = [self.get_contact_by_id(p,short_id,cat=sc) for p in other_prgs]
        exp_res = None
        exp_fuzzy = True
        for x in other_results:
            c = len([y for y in other_results if y==x])
            if c>=2 or (sc in ['base-ribose','base-phosphate'] and c>=1):
                assert x==exp_res or (exp_res is None)
                exp_fuzzy = False
                exp_res = x
        if exp_res is None:
            exp_res = ""
            
        res = None
        if prg_res != "" and prg_res == exp_res:
            res = "tp"
        elif prg_res == "" and prg_res == exp_res:
            res = "tn"
        elif prg_res != "" and prg_res != exp_res:
            res = "fp"
        elif prg_res == "" and prg_res != exp_res:
            res = "fn"

        assert res in ('tp','tn','fp','fn')
            
        if details:
            assert prg in ['CL']
            assert self.dd is not None
            n_type = self.dd.get_n_type(d_id)
            res_sc = ""
            comment = ""
            if res=='fp':
                res_sc,comment = self.classify_prg_result_details_fp(short_id,sc,prg_res,exp_res,other_results)
            return {"res":res,"res_sc":res_sc,"d_id":d_id,"n_type":n_type,
                    "cl_res":prg_res,"exp_res":exp_res,"exp_fuzzy":exp_fuzzy,
                    "other_results":zip(other_prgs,other_results),
                    "extra":comment}

        return res

def combine_graphs_new_eval(options):
    filter_doublets = None
    if options.only_doublets_from:
        filter_doublets = load_doublets_from_file(options)

    verbosity = int(options.verbosity)

    CATEGORIES = ["bp-classic","bp-non-classic","stacking","base-phosphate","base-ribose"]
    PRG_BY_CAT = {
        "bp-classic": ["CL","RV","MC","FR"],
        "bp-non-classic": ["CL","RV","MC","FR"],
        "stacking": ["CL","MC","MO","FR"],
        "base-phosphate": ["CL","FR"],
        "base-ribose": ["CL","FR"],
    }
    STEPS = np.linspace(0.0, 1.0, 21)

    dd = DoubletsDict(options.data_dir,reduced_atoms='*')
    dd.load_pdb(options.pdb_id)
    
    graphs = GraphsCollection(
                pdb_id=options.pdb_id,
                cl_fn=options.classifier_graph,
                fr_fn=options.fr3d_graph,
                mc_fn=options.mc_annotate_graph,
                mo_fn=options.moderna_graph,
                rv_fn=options.rnaview_graph,
                dd=dd,
                filter_doublets=filter_doublets)
    groups = GroupsTool(options.groups,filter_doublets=filter_doublets)
    
    res = {}
    for sc in CATEGORIES:
        P_LIST = PRG_BY_CAT[sc]
        
        for prg in sorted(P_LIST):
            OTHER_PRGS = [x for x in P_LIST if x!=prg and x!='CL']
            if len(OTHER_PRGS)==0:
                continue
            all_ids = set(groups.get_all_ids())
            if sc not in ['bp-classic','bp-non-classic','stacking']:
                all_ids = all_ids.union([GroupsTool.reverse_id(x) for x in groups.get_all_ids()])
                            
            for d_id in all_ids:
                keys = []
                keys_d = []
            
                short_id = re.sub("^[^:]*:", "", d_id)
                k = graphs.classify_prg_result(prg,short_id,sc)
                keys.append("evaluation/%(prg)s/%(k)s/%(sc)s" % locals())
                if prg in ['CL']:
                    k_details = graphs.classify_prg_result(prg,short_id,sc,details=True)
                    cl_res = k_details['cl_res']
                    exp_res = k_details['exp_res']
                    n_type = k_details['n_type']
                    sc_simple = sc.replace("bp-classic","bp").replace("bp-non-classic","bp")
                    if not k_details['exp_fuzzy']:
                        keys.append("evaluation/%(prg)s_ignore_fuzzy/%(k)s/%(sc)s" % locals())
                    for i in STEPS:
                        kk = graphs.classify_prg_result(prg,short_id,sc,min_contact_weight=i)
                        keys.append("evaluation/%(prg)s/%(kk)s/%(sc)s/%(i).2f" % locals())
                    
                    if k in ['tp','fp']:
                        keys.append("evaluation/%(prg)s/%(k)s/%(sc_simple)s/%(cl_res)s/%(n_type)s" % locals())
                        if k=='fp':
                            kk = k+"-"+k_details['res_sc']
                            keys.append("evaluation/%(prg)s/%(kk)s/%(sc_simple)s/%(cl_res)s/%(n_type)s" % locals())
                            keys_d.append(("evaluation-doublets/%(prg)s/%(kk)s/%(sc_simple)s/%(cl_res)s/%(n_type)s" % locals(),d_id))
                            keys_d.append(("evaluation-details/%(prg)s/%(kk)s/%(sc_simple)s" % locals(),k_details))
                        else:
                            keys_d.append(("evaluation-doublets/%(prg)s/%(k)s/%(sc_simple)s/%(cl_res)s/%(n_type)s" % locals(),d_id))
                    elif k=='fn':
                        keys.append("evaluation/%(prg)s/%(k)s/%(sc_simple)s/%(exp_res)s/%(n_type)s" % locals())
                        keys_d.append(("evaluation-doublets/%(prg)s/%(k)s/%(sc_simple)s/%(exp_res)s/%(n_type)s" % locals(),d_id))

                for key in keys:
                    if not res.has_key(key):
                        res[key] = 0
                    res[key] += 1
                for key,d in keys_d:
                    if not res.has_key(key):
                        res[key] = []
                    res[key].append(d)

        if verbosity>0:
            for i in STEPS:
                row = {}
                for key in ['tp','fp','tn','fn']:
                    full_key = 'evaluation/CL/%s/%s/%.2f'%(key,sc,i)
                    row[key] = res.get(full_key,0)
                
                row = confusion_matrix_params(row)
                if i==0.5:
                    print "==================================================="
                print "sc=%s[%.2f]:" % (sc,i),
                print "tp=%(tp)d tn=%(tn)d fp=%(fp)d fn=%(fn)d tpr=%(tpr).3f fpr=%(fpr).3f acc=%(acc).3f spc=%(spc).3f mcc=%(mcc).3f" % row
                if i==0.5:
                    print "==================================================="


    save_json(options.output_json, res)


def main():
    (parser, options, args) = parse_args()
    if options.new_eval_mode:
        assert options.groups is not None
        assert options.rnaview_graph is not None
        assert options.mc_annotate_graph is not None
        assert options.moderna_graph is not None
        assert options.fr3d_graph is not None
        assert options.classifier_graph is not None
        assert options.pdb_id is not None
        combine_graphs_new_eval(options)
    elif options.correlation:
        assert options.rnaview_graph is not None
        assert options.mc_annotate_graph is not None
        assert options.moderna_graph is not None
        assert options.fr3d_graph is not None
        assert options.classifier_graph is not None
        assert options.pdb_id is not None
        combine_graphs_correlation(options)
    elif options.groups:
        assert options.rnaview_graph is None
        assert options.mc_annotate_graph is None
        assert options.moderna_graph is None
        assert options.fr3d_graph is None
        assert options.pdb_id is not None
        combine_graphs_using_groups(options)
    elif options.eval_mode:
        assert options.rnaview_graph is not None
        assert options.mc_annotate_graph is not None
        assert options.fr3d_graph is not None
        assert options.classifier_graph is not None
        assert options.pdb_id is not None
        combine_graphs_eval(options)
    else:
        combine_graphs_old(options)
    
if __name__ == '__main__':
    main()

