#!/usr/bin/env python
import os
from utils import *
from distances import rmsd_distance, bp_distance, bph_distance, doublet_params_dict


R_ATOMS = ['N1','N2','N3','N4','N6','N7','P']
R_ATOMS += ['C2','C4','C5','C6','C8',"C1'"]
R_ATOMS += ['O2']
R_ATOMS += ["O2'","O3'","O4'"]
R_ATOMS += ["OP1","OP2","O5'","NEXT:O3'"]
R_ATOMS += ['N1','C6','O6','C5','C4','N3','C2','N2','N7','C8','N9']

def compute_doublet_params(ids, params_type):
    dp = __import__("doublet-params")

    dd = DoubletsDict("gc-data/", reduced_atoms=R_ATOMS)
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
            print "INVALID doublet! %s" % (d_id)
            continue
        res[d_id] = d
        pbar.update(i)
    pbar.finish()
    return res

def filter_by_params(doublets, limits, all_params_dict):
    res = []
    for d_id in doublets:
        if not all_params_dict.has_key(d_id):
            continue
        p = all_params_dict[d_id]
        ok = all([p[k]>=limit['min'] and p[k]<=limit['max'] for k,limit in limits.items()])
        if ok:
            res.append(d_id)
    return res

def load_unclassified_groups(n_type):
    tmp_json = "/tmp/unclassified.json"
    os.system("wget -q 'http://rnacontacts-vm/unclassified-groups/json?n_type=%s&search=' -O '%s'" % (n_type,tmp_json))
    return load_json(tmp_json)

def download_members(id):
    url = "http://rnacontacts-vm/doublets/group/info/%d/download-members/" % id
    tmp_json = "/tmp/members.json"
    os.system("wget -q '%s' -O '%s'" % (url,tmp_json))
    res = load_json(tmp_json)
    os.unlink(tmp_json)
    return res

def normalize_members(members,type='bp'):
    sup_atoms = (SELECTED_ATOMS,[])
    rmsd_atoms = ([],SELECTED_ATOMS)
    if type=='base-ribose':
        rmsd_atoms = ([], ["O2'","O3'","O4'","C1'","C2'","C3'","C4'"])
    elif type=='ribose-ribose':
        sup_atoms = (["O2'","O3'","O4'","C1'","C2'","C3'","C4'"],[])
        rmsd_atoms = ([], ["O2'","O3'","O4'","C1'","C2'","C3'","C4'"])
    
    
    dd = DoubletsDict("gc-data",reduced_atoms=sum(sup_atoms,[])+sum(rmsd_atoms,[]))
    dd.load_pdb_files(members,verbose=True)
    ref = dd.get(members[0])
    new_members = []
    for d_id in members:
        (p1,p2) = dd.get(d_id)
        dist1 = rmsd_distance((p1,p2),ref,(SELECTED_ATOMS,[]),([],SELECTED_ATOMS))
        dist2 = rmsd_distance((p2,p1),ref,(SELECTED_ATOMS,[]),([],SELECTED_ATOMS))
        if dist1>dist2:
            print "reversing: %s" %d_id
            d_id = GroupsTool.reverse_id(d_id)
        new_members.append(d_id)
    return new_members

def add_group(name,elems,groups,type="bp"):
    selected_groups = [groups[x-1] for x in elems]
    print name, set(sum([x['notes'].split(",") for x in selected_groups],[])), [x['id'] for x in selected_groups]
    for i,x in enumerate(selected_groups,start=1):
        members = download_members(x['id'])
        members = normalize_members(members,type)
        out_json = "/tmp/%s%d.json"%(name,i)
        out_pdb = "/tmp/%s%d.pdb.gz"%(name,i)
        save_json(out_json,members)
        os.system("./gen-pdb.py --type='%s' --input-json='%s' -o '%s'" % (type,out_json,out_pdb))
        if type=='base-ribose':
            out_pdb = re.sub(".pdb(.gz)$",r"-centers.pdb\1",out_pdb)
            os.system("./gen-pdb.py --type='ribose-centers' --input-json='%s' -o '%s'" % (out_json,out_pdb))
        

def process_pu_pu_groups(groups):
    add_group("base-ribose-stacking",[25,50,86,98],groups,type="base-ribose")
    add_group("ribose-ribose",[64,65,72,73,74,79],groups,type="ribose-ribose")
    add_group("diagonal a",[1,4,5,6,7,9],groups)
    add_group("diagonal b",[60,61,63,69,82,84,87,88],groups)
    add_group("intercal",[2,13,16],groups)
    add_group("diagonal c",[3],groups)

def find_diagonal_groups():
    result = {}
    for pdb_id in read_file("data-all.txt").strip().split("\n"):
        gr_fn = PDBObject(pdb_id,"contacts_MC")
        graph = GraphTool(gr_fn)
        nodes_with_data = sorted(list(graph.nodes(data=True)), key=lambda x: ((x[0][0],int(x[0][1:])),x[1]))
        nodes = [x[0] for x in nodes_with_data]
        nodes_data = [x[1] for x in nodes_with_data]
        n_types = dict([(nodes[i],nodes_data[i]['resname']) for i in xrange(len(nodes))])
        
        
        def add(label,r1,r2):
            n_type = n_types[r1]+n_types[r2]
            if not result.has_key(label):
                result[label] = {}
            if not result[label].has_key(n_type):
                result[label][n_type]=set()
            full_id = "%s:%s:%s" % (pdb_id.upper(),r1,r2)
            result[label][n_type].add(full_id)
            print "adding %s to %s/%s" % (full_id,label,n_type)
        
        for i,rid1 in enumerate(nodes):
            rid2_list = []
            if i>0:
                rid2_list.append(nodes[i-1])
            if i<len(nodes)-1:
                rid2_list.append(nodes[i+1])
            for rid2 in rid2_list:
                if rid1[0]==rid2[0] and int(rid1[1:])==int(rid2[1:])-1:
                    for rid3 in graph.g.neighbors(rid1):
                        if rid3==rid1 or rid3==rid2:
                            continue
                        
                        if graph.get_contact_by_id(rid1+":"+rid3,cat="stacking")!='' and \
                           graph.get_contact_by_id(rid2+":"+rid3,cat="bp")!='' and \
                           graph.get_contact_by_id(rid1+":"+rid2,cat="stacking")=='':
                            add("diagonal-c",rid1,rid2)
                        if graph.get_contact_by_id(rid1+":"+rid3,cat="stacking")!='' and \
                           graph.get_contact_by_id(rid2+":"+rid3,cat="stacking")!='':
                            add("long-stacking-c",rid1,rid2)
                        if graph.get_contact_by_id(rid1+":"+rid3,cat="bp")!='':
                            for rid4 in graph.g.neighbors(rid2):
                                if rid4==rid1 or rid4==rid2 or rid4==rid3:
                                    continue
                                if graph.get_contact_by_id(rid3+":"+rid4,cat="stacking")!='' and \
                                   graph.get_contact_by_id(rid4+":"+rid2,cat="stacking")!='':
                                    add("long-diagonal-c",rid1,rid2)
                        if graph.get_contact_by_id(rid1+":"+rid2,cat="stacking")!='' and \
                           graph.get_contact_by_id(rid1+":"+rid3,cat="bp") in ['WW_cis']:
                            for rid4 in graph.g.neighbors(rid2):
                                if rid4==rid1 or rid4==rid2 or rid4==rid3:
                                    continue
                                if graph.get_contact_by_id(rid3+":"+rid4,cat="stacking")!='' and \
                                   graph.get_contact_by_id(rid2+":"+rid4,cat="bp") in ['WW_cis'] and \
                                   graph.get_contact_by_id(rid1+":"+rid4,cat="bp")=='' and \
                                   graph.get_contact_by_id(rid2+":"+rid3,cat="bp")=='' and \
                                   graph.get_contact_by_id(rid1+":"+rid4,cat="stacking")=='':
                                    add("diagonal-nc-ww",rid1,rid4)
    for label,v in result.items():
        limit_consecutive = None
        limits = None
        if label=='long-stacking-c':
            limit_consecutive = True
            limits = {'nn_ang_norm': {'min':0.0,'max':40.0}, 
                      'dist': {'min':6.4,'max':9.0}, 
                      'n12cc_ang': {'min': 0.0,'max':45.0}, 
                      }
        elif label=='diagonal-c':
            limit_consecutive = True
            limits = {'nn_ang_norm': {'min':0.0,'max':60.0}, 
                      'dist': {'min':4.0,'max':9.0}, 
                      'n12cc_ang': {'min': 30.0,'max':70.0}, 
                      }
        elif label=='diagonal-nc-ww':
            limit_consecutive = False
            limits = {'nn_ang_norm': {'min':0.0,'max':40.0}, 
                      'dist': {'min':5.5,'max':10.5}, 
                      'dist_z': {'min':1.8,'max':5.5}, 
                      'n12cc_ang': {'min': 50.0,'max':80.0}, 
                      }
        if label not in ['diagonal-nc-ww','diagonal-c','long-stacking-c']:
            continue
            
        for n_type,members in v.items():
            print label,n_type,len(members)

            if limits is not None:
                all_params = compute_doublet_params(members, "bp")
                members = filter_by_params(members, limits, all_params)
                print "after filtering: %d" % len(members)

            out_json = "expert/ref/other_%s_%s.json"%(n_type,label)
            out_pdb = "expert/ref/other_%s_%s.pdb.gz"%(n_type,label)
            save_json(out_json,sorted(list(members)))
            os.system("./gen-pdb.py --type='%s' --input-json='%s' -o '%s' --limit=500" % (type,out_json,out_pdb))


# process_pu_pu_groups(load_unclassified_groups("Pu-Pu"))
find_diagonal_groups()
