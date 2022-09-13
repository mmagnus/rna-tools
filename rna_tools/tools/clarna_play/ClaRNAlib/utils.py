#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
import os
import re
import math
import simplejson as json
import urllib.request, urllib.parse, urllib.error
import io
import gzip
import pickle
import hashlib
import itertools
from datetime import datetime
from optparse import OptionParser

from numpy import array
from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial import KDTree
import tempfile

from rna_tools.tools.clarna_play.ClaRNAlib.distances import normalize_points, fit_points
from rna_tools.tools.clarna_play.ClaRNAlib.cl_settings_global import *

# deprecated
FR3D_URL = "http://rna.bgsu.edu/FR3D/AnalyzedStructures"
FR3D_CACHE = "/home/twalen/fr3d"

N_TYPES = [n1+n2 for n1,n2 in itertools.product(["A","C","U","G"],repeat=2)]

CL_CAT = {
    "bp":[
        "HH_cis",
        "HH_tran",
        "HS_cis",
        "HS_tran",
        "HW_cis",
        "HW_tran",
        "SH_cis",
        "SH_tran",
        "SS_cis",
        "SS_tran",
        "SW_cis",
        "SW_tran",
        "WH_cis",
        "WH_tran",
        "WS_cis",
        "WS_tran",
        "WW_cis",
        "WW_tran"
    ],
    "stacking":[
        "<<",">>","<>","><"
    ],
    "base-phosphate":[
        "H_0BPh",
        "SW_2BPh",
        "S_1BPh",
        "W_345BPh",
        "W_6BPh",
        "H_789BPh",
    ],
    "base-ribose":[
        "H_0BR",
        "S_1BR",
        "SW_2BR",
        "W_345BR",
        "W_6BR",
        "H_789BR",
    ],
    "other":[
        "diagonal-c",
        "diagonal-nc-ww",
        "long-stacking-c",
    ],
    "other2":[
        "base-ribose-stacking",
    ],
    "other3":[
        "phosphate-stacking1",
        "phosphate-stacking2",
    ],
}


SELECTED_ATOMS = ["C1'","C4","C6","C2"]

# minimal list of atoms required to be present in VALID residues
REQ_ATOMS_LIST = ('C2','C4','C6',"C1'",'P',"O3'")

def status_msg(msg,new_line=True):
    if new_line:
        print(msg)
    else:
        print(msg, end=' ')
    sys.stdout.flush()

def bench_start(msg,new_line=False):
    global BENCH_T
    BENCH_T = datetime.now()
    status_msg(msg + " .. ",new_line=new_line)

def bench_stop(msg=""):
    t = datetime.now()
    diff = t-BENCH_T
    if msg!="":
        prefix = msg + " .. "
    else:
        prefix = ""
    diff_sec = float(diff.seconds)+float(diff.microseconds)/1000000
    status_msg("%sDONE [%.3f sec.]"%(prefix,diff_sec))

def write_file(fn, contents):
    if re.match(".*.gz$",fn):
        f = gzip.open(fn,"w")
    else:
        f = open(fn, "w")
    f.write(contents)
    f.close()

def find_files_in_dir(dirname):
    """find all files in directory dirname (and its subdirectories)"""
    res = []
    for root, _dirs, files in os.walk(dirname):
        for f in files:
            res.append(os.path.join(root, f))
    return res

def ensure_dir(f):
    # from http://stackoverflow.com/questions/273192/python-best-way-to-create-directory-if-it-doesnt-exist-for-file-write
    d = os.path.dirname(f)
    if d!='' and not os.path.exists(d):
        os.makedirs(d)

def pdb_files(dirname,include_positive=True,include_negative=False):
    res = []
    for f in find_files_in_dir(dirname):
        id = os.path.basename(f)
        if re.match(".*merged.*", id):
            continue
        if re.match(".*negative.*", id):
            if include_negative:
                res.append(f)
        elif re.match(r"^.*\.pdb$", id):
            if include_positive:
                res.append(f)
    return res

def simplify_residue(r):
    res = {}
    for a in r:
        if a.id in res:
            print("duplicate atom a.id=%s" % (a.id))
            return None
        assert((a.id in res) == False)
        res[a.id] = a.get_coord()
    return res

def get_atoms(structure):
    """NOT USED"""
    SELECTED_ATOMS = ["C1'","N1","N3"]
    result = []
    for r in structure.get_residues():
        res = []
        for atom_name in SELECTED_ATOMS:
            aa = None
            for a in r:
                if a.id==atom_name:
                    aa = a.get_coord()
                    break
            if aa is None:
                res = None
                break
            else:
                res.append(aa)
        result.append(res)
    if len(result)==0:
        return None
    return array(result,'f')


def simplify_structure(s):
    res = {}
    if hasattr(s, 'get_chains'):
        l = s.get_chains()
    else:
        l = s
    for (chain_num,chain) in enumerate(l):
        for r in chain:
            rr = simplify_residue(r)
            if "C1'" in rr:
                cat = 0
            elif "O3'" in rr:
                cat = 1
            else:
                cat = 2
            for a in list(rr.keys()):
                key = "%d%d-%s" % (chain_num,cat, a)
                res[key]=rr[a]
    return res

def load_clusters(fn):
    f = open(fn)
    res = json.load(f)
    f.close()
    return res

def parse_pdb(s):
    f = io.StringIO(s)
    parser = PDB.PDBParser()
    return parser.get_structure("c",f)

def load_pdb(fn):
    if re.match(".*.gz$",fn):
        f = gzip.open(fn)
    else:
        f = open(fn)
    parser = PDB.PDBParser()
    res = parser.get_structure("c",f)
    for a in res.get_atoms():
        if re.match(r'^[A-Z]{1,2}[0-9]?\*$',a.id):
            a.id = a.id.replace("*","'")
    return res

def save_pdb(fn,structure):
    ensure_dir(fn)
    if re.match(".*.gz$",fn):
        f = gzip.open(fn,"w")
    else:
        f = open(fn,"w")
    f.write(structure2string(structure))
    f.close()

def structure2string(structure):
    f = StringIO()
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(f)
    return f.getvalue()

def string2structure(pdb_str):
    f = StringIO(pdb_str)
    parser = PDB.PDBParser()
    return parser.get_structure("structure", f)

def load_json(fn):
    if re.match(".*.gz$",fn):
        f = gzip.open(fn)
    else:
        f = open(fn)
    return json.load(f)

def save_json(fn,data,indent=None):
    ensure_dir(fn)
    if re.match(".*.gz$",fn):
        f = gzip.open(fn,"w")
    else:
        f = open(fn,"w")
    return json.dump(data,f,indent=indent)

def cache_init():
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)

def cache(prefix):
    """cache decorator"""
    def add_cache(f):
        def new_f(*args):
            cache_key = hashlib.md5(str(args)).hexdigest()
            file_name = os.path.join(CACHE_DIR, prefix + '_' + cache_key)

            try:
                return pickle.load(open(file_name))
            except IOError:
                value = f(*args)
                pickle.dump(value, open(file_name, 'w'))
            return value
        return new_f

    return add_cache

def pickle_save(fn, value):
    ensure_dir(fn)
    pickle.dump(value, open(fn, 'w'))

def pickle_load(fn):
    if not os.path.isfile(fn):
        return None
    return pickle.load(open(fn))

@cache("get_pdb")
def _get_pdb_web_page(pdb_id):
    pdb_id = pdb_id.lower()
    assert len(pdb_id)==4
    f = urllib.request.urlopen("http://www.pdb.org/pdb/explore/explore.do?structureId="+pdb_id.lower())
    data = f.read()
    f.close()
    return data

def get_pdb_info(pdb_id):
    pdb_id = pdb_id.lower()
    res = {"pdb_id":pdb_id}
    if re.match(r"^zz",pdb_id):
        return res
    data = _get_pdb_web_page(pdb_id)
    for field,regular_expr in [
        ('release_date',r'<span class="se_key">Release:.*<td class="perc33">\s+([^\s]+)\s+</td>'),
        ('title',r'<div id="se_structureTitle".*<span class="h3">\s*([^<]*)\s*</span>'),
        ('resolution',r'<td id="se_xrayResolution">\s*([^\s+]*)\s*</td>'),
        ('em_resolution',r'<span class="se_key">EM Resolution.*</span>\s*([\d\.]+)\s*')
        ]:
        regex = re.compile(regular_expr, re.MULTILINE|re.DOTALL)
        m = regex.search(data)
        if m:
            res[field] = m.group(1)
    if 'em_resolution' in res and 'resolution' not in res:
        res['resolution'] = res['em_resolution']
    return res

def load_graph(filename):
    from networkx.readwrite import json_graph
    if re.match(r".*\.gz$",filename):
        f = gzip.open(filename,"r")
    else:
        f = open(filename,"r")
    return json_graph.load(f)

def save_graph(filename,g,indent=None):
    from networkx.readwrite import json_graph
    print('json save to file:', filename)
    if re.match(r".*\.gz$",filename):
        f = gzip.open(filename,"w")
    else:
        f = open(filename,"w")
    f.write(json.dumps(json_graph.node_link_data(g)))
    j = json_graph.node_link_data(g)
    t = ''
    for i in j['links']:
        # ?<> 
        if i['type'] != 'dist':
            if '?' not in i['key']:
                # print(i)
                if i['key'].startswith('WW') or '>' in i['key'] or '<' in i['key']:
                   t += i['source'].ljust(10) + i['target'].ljust(10) + i['key'].ljust(10) + str(round(i['weight'], 2)) + '\n'
    print(t)
    return json_graph.node_link_data(g)#, f, indent=indent)


def load_motifs(dirname,include_positive=True,include_negative=False):
    parser = PDB.PDBParser()
    res = {}
    for f in pdb_files(dirname,include_positive,include_negative):
        ff = os.path.basename(f)
        key1 = ff[0:2]
        key2 = re.sub(r"^(..)_(desc|negative)_(.*)_from.*$",r"\3",ff)
        if key1 not in res:
            res[key1] = {}
        if key2 not in res[key1]:
            res[key1][key2]=[]
        # print "loading %s/%s" % (key1,key2)
        s = parser.get_structure("c",f)
        s = simplify_structure(s,key1)
        if s is None:
            print("ignoring: %s" % f)
            continue
        res[key1][key2].append(s)
    return res

def dist(p1, p2):
    return math.sqrt(sq_dist(p1, p2))

def sq_dist(p1, p2):
    assert(len(p1)==len(p2))
    assert(len(p1)==3)
    return ((p1[0]-p2[0])*(p1[0]-p2[0]) +
            (p1[1]-p2[1])*(p1[1]-p2[1]) +
            (p1[2]-p2[2])*(p1[2]-p2[2]))

def res_p_distance(r1, r2):
    if 'P' not in r1 or 'P' not in r2:
        return 1000
    return dist(r1['P'], r2['P'])

def res_atom_distance(r1, r2, from_atom, to_atom=None):
    if to_atom is None:
        to_atom = from_atom
    if from_atom not in r1 or to_atom not in r2:
        return 1000
    return dist(r1[from_atom], r2[to_atom])


def res_c1p_distance(r1, r2):
    if "C1'" not in r1 or "C1'" not in r2:
        return 1000
    return dist(r1["C1'"], r2["C1'"])

def res_c2_distance(r1, r2):
    if "C2" not in r1 or "C2" not in r2:
        return 1000
    return dist(r1["C2'"], r2["C2'"])

def read_file(fn):
    f = open(fn)
    res = f.read()
    f.close()
    return res

def save_file(fn, contents):
    ensure_dir(fn)
    f = open(fn, "w")
    f.write(contents)
    f.close()

def run_rnaview(fn):
    os.putenv("RNAVIEW", RNAVIEW_HOME)
    os.system("%s/bin/rnaview %s > /dev/null 2> /dev/null" % (RNAVIEW_HOME,fn))
    out_fn = fn + ".out"
    res = read_file(out_fn)
    os.unlink(out_fn)
    return res

def run_mc_annotate(fn):
    out_fn = fn + ".mc_out"
    os.system("%s/MC-Annotate %s > %s" % (MCA_HOME,fn,out_fn))
    res = read_file(out_fn)
    os.unlink(out_fn)
    return res

def run_moderna(fn):
    out_fn = fn + ".mc_out"
    os.system("%s %s/moderna-stackings.py --use-fr3d-format %s > %s" % (PYTHON_CMD, MODERNA_STACKINGS_HOME, fn,out_fn))
    res = read_file(out_fn)
    os.unlink(out_fn)
    return res


def run_new_alg2(fn,db="rnaview"):
    out_fn = fn + ".my_out"
    os.system("./find-contacts.py --new-alg2 --input=%s --db=contacts-db/%s.setA.merged > %s" % (fn,db,out_fn))
    res = read_file(out_fn)
    os.unlink(out_fn)
    return res

def run_new_alg3(fn,db="rna-view.json"):
    out_fn = fn + ".my_out"
    os.system("%s ./clarna.py --input=%s --lib=%s > %s" % (PYTHON_CMD,fn,db,out_fn))
    res = read_file(out_fn)
    os.unlink(out_fn)
    return res

def run_fr3d_old(fn):
    pdb_id = os.path.basename(fn)[0:4].upper()
    result = ""
    for t in ['basepairs','stacking','base_phosphate']:
        fn = os.path.join(FR3D_CACHE,pdb_id + "_" + t + ".html")
        if os.path.isfile(fn):
            print("using cache")
            f = open(fn)
        else:
            url = FR3D_URL + "/" + pdb_id + "/" + pdb_id + "_" + t + ".html"
            f = urllib.request.urlopen(url)
        html = f.read()
        f.close()
        # extract pre fragment
        result += "# %s\n" % t
        in_pre = False
        for line in html.split("\n"):
            if in_pre:
                if re.match(r"^\s*</pre>",line):
                    in_pre = False
                else:
                    m = re.match(r"^.*\((.)\).*\((.)\)",line)
                    if m:
                        result += line + "\n"
            else:
                if re.match(r"^\s*<pre>",line):
                    in_pre = True
    return result

def run_fr3d(fn):
    fr3d_home = FR3D_HOME
    tmp_fn = None
    if re.match("^.*.gz$",fn):
        (_fd,tmp_fn) = tempfile.mkstemp(".fr.pdb")
        retcode = os.system("cat '%(fn)s' | gunzip -c > '%(tmp_fn)s'" % locals())
        assert retcode==0, "error on executing gunzip"
        pdb_fn = tmp_fn
    else:
        pdb_fn = fn
    (_fd2,out_fn) = tempfile.mkstemp(".fr.csv")
    cmd = """octave --eval "cd '%(fr3d_home)s'; addpath(genpath('%(fr3d_home)s')); annotate_pdb('%(pdb_fn)s','%(out_fn)s')" """ % locals()
    retcode = os.system(cmd)
    assert retcode==0, "error on executing fr3d!"
    res = read_file(out_fn)
    os.unlink(out_fn)
    if tmp_fn is not None:
        os.unlink(tmp_fn)
    return res

def download_rna_hub(pdb_id):
    cache_fn = os.path.join(FR3D_CACHE,pdb_id.lower() +".csv")
    if os.path.isfile(cache_fn):
        print("using cache")
        f = open(cache_fn)
    else:
        url = "http://rna.bgsu.edu/rna3dhub/pdb/%s/interactions/fr3d/all/csv" % pdb_id.upper()
        print(url)
        f = urllib.request.urlopen(url)
    res = f.read()
    f.close()
    new = []
    for line in res.split("\n"):
        if line=="":
            continue
        if "No interactions" in line:
            break
        tmp = line.replace('"',"").replace("__","").split(",")
        assert len(tmp)==6
        (_pdb1,model1,chain1,r1,num1) = tmp[0].split("_")[0:5]
        (_pdb2,model2,chain2,r2,num2) = tmp[5].split("_")[0:5]
        name1 = "%s_%s_%s_%s" % (model1,chain1,num1,r1)
        name2 = "%s_%s_%s_%s" % (model2,chain2,num2,r2)
        for i in (1,2,3,4):
            if tmp[i]!='':
                new.append('"%s","%s","%s"' % (name1,tmp[i],name2))
    return "\n".join(new)

def simplify_desc(desc):
    desc = re.sub(r"syn"," ",desc)
    desc  = re.sub(r"[/!()\?]","",desc)
    desc = desc.strip()
    desc = re.sub(r"\s+","_",desc)
    desc  = desc.replace("+","P")
    desc  = desc.replace("-","M")
    desc = desc.replace("stacked", "stacking")
    return desc

def simplify_desc_mc(desc):
    desc = desc.strip()
    desc = re.sub(r"\s+","_",desc)
    desc = desc.replace("+","P")
    desc = desc.replace("-","M")
    desc = desc.replace(">","gt")
    desc = desc.replace("<","lt")
    desc = re.sub(r"[/!()\?]","",desc)
    return desc

def mc_to_rv_desc(desc):
    if desc is None or desc=='':
        return None
    if re.match('^stacking',desc):
        return desc
    desc = desc.replace("O2P","H")
    desc = desc.replace("O2'","S")
    desc = desc.replace("Bs","S")
    desc = desc.replace("Bh","W")
    m = re.match(r"^([sHSW]).([sHSW]).*_(cis|tran).*$",desc)
    if m:
        return m.group(1)+m.group(2)+"_"+m.group(3)
    else:
        return "UNK(%s)" % desc

def parse_mc_annotate(output,base_pairs=True,stackings=False,use_chain=False):
    res = []
    residues = {}
    in_residue_conformations = False
    in_bp = False
    in_stackings = False
    for line in output.split("\n"):
        if in_residue_conformations:
            # C956 : A C3p_endo anti
            # C957 : A C3p_endo syn
            m = re.match("^'?(.)'?(\d+) : (.) ([^ ]*) ([^ ]*)$", line)
            if m:
                chain1 = m.group(1)
                num1 = m.group(2)
                n_type = m.group(3)
                conf = m.group(5)
                if use_chain:
                    num1 = chain1+num1
                residues[num1] = {'resname': n_type, 'conf': conf}
        if in_bp and base_pairs:
            # C902-C970 : G-C Ww/Ww pairing antiparallel cis one_hbond 130
            # '9'77-'9'102 : G-A Ww/O2P pairing
            m = re.match("^'?(.)'?(\d+)-'?(.)'?(\d+) : (.)-(.) (.*)$", line)
            if m:
                chain1 = m.group(1)
                num1 = m.group(2)
                chain2 = m.group(3)
                num2 = m.group(4)
                if use_chain:
                    num1 = chain1+num1
                    num2 = chain2+num2
                r1 = m.group(5)
                r2 = m.group(6)
                desc = m.group(7).strip()
                desc2 = simplify_desc_mc(desc)
                # print in_stackings,stackings,num1,num2,r1,r2,desc2
                res.append([(num1,num2),(r1,r2),desc2])
        if in_stackings and stackings:
            # C968-C969 : adjacent_5p upward
            m = re.match("^'?(.)'?(\d+)-'?(.)'?(\d+) : (.*)$", line)
            if m:
                chain1 = m.group(1)
                num1 = m.group(2)
                chain2 = m.group(3)
                num2 = m.group(4)
                if use_chain:
                    num1 = chain1+num1
                    num2 = chain2+num2
                r1 = '-'
                r2 = '-'
                desc = m.group(5).strip()
                desc2 = "stacking_"+simplify_desc_mc(desc)
                res.append([(num1,num2),(r1,r2),desc2])
        if (re.match("^(Non-)?Adjacent stackings ---.*$",line) or
            re.match("^Stackings ---.*$",line)): # for moderna-stackings
            in_stackings = True
            in_bp = False
            in_residue_conformations = False
        elif re.match("^Base-pairs ---.*$",line):
            in_bp = True
            in_stackings = False
            in_residue_conformations = False
        elif re.match("^Residue conformations ---.*$",line):
            in_residue_conformations = True
            in_bp = False; in_stackings = False
    return res, residues

def parse_rnaview(output, my_output=False, base_pairs=True,stackings=False,use_chain=False,cut_desc=True):
    # sample
    #           1         2         3         4
    # 012345678901234567890123456789012345678901234567890123456789
    #      2_70, C:     1 G-C    69 C: +/+ cis         XIX
    res = []
    if my_output:
        in_bp = True
    else:
        in_bp = False
    for line in output.split("\n"):
        if "!(" in line or "!1" in line:
            continue
        if in_bp and len(line)>33:
            if not stackings and "stacked" in line:
                continue
            if not base_pairs and not "stacked" in line:
                continue
            c1 = line[11]
            c2 = line[30]
            num1 = line[13:19].strip()
            num2 = line[23:29].strip()
            if use_chain:
                num1 = c1+num1
                num2 = c2+num2
            r1 = line[20];
            r2 = line[22]
            if my_output:
                if cut_desc:
                    desc = line[33:54].strip()
                else:
                    desc = line[33:].strip()
                desc2 = desc
            else:
                desc = line[33:].strip()
                desc2 = simplify_desc(desc)
            # print "num1=#%s# num2=#%s# r=#%s#%s# desc=#%s/%s#" % (num1,num2,r1,r2,desc,desc2)
            res.append([(num1,num2),(r1,r2),desc2])
        if re.match("BEGIN_base-pair", line):
            in_bp = True
        elif re.match("END_base-pair", line):
            in_bp = False
    return res

def parse_fr3d_old(output,base_pairs=True,stackings=False,use_chain=True):
    in_bp = False
    in_stackings = False
    res = []
    for line in output.split("\n"):
        if re.match(r'^# base',line):
            in_bp = True
            in_stackings = False
        elif re.match(r'# stacking',line):
            in_stackings = True
            in_bp = False
        if len(line)>20:
            m = re.match(r"^\s*\d+\s+(.)\s*([\dA-Z]+)\((.)\)\s+-\s+(.)\s*([\dA-Z]+)\((.)\)\s+-\s*([^\s]+)\s*-",line)
            if m is None:
                continue
            #           1         2         3         4
            # 0123456789012345678901234567890123456789012
            # ------------------------------------------
            #     1 C  14(A) - U  13(A) -  H_9BPh -    0
            # or
            #    1 G   9(A) - C  25(A) -    cWW  -    0
            # or
            #   2210 G1003A(A) - G1003(A) -    s53  -    0
            r1 = m.group(1)
            num1 = m.group(2).strip()
            chain1 = m.group(3)
            r2 = m.group(4)
            num2 = m.group(5).strip()
            chain2 = m.group(6)
            desc = m.group(7).strip()
            if in_stackings:
                desc = "stacking_"+desc
            name1 = num1
            name2 = num2
            if use_chain:
                name1 = chain1+num1
                name2 = chain2+num2
            else:
                name1 = num1
                name2 = num2

            if (stackings and in_stackings) or (base_pairs and in_bp):
                res.append([(name1,name2),(r1,r2),desc])
    return res

def parse_fr3d(output,base_pairs=True,stackings=False,use_chain=True):
    res = []
    for line in output.strip().split("\n"):
        # example:
        # "1_B_1_G","s53","1_B_0_G"
        tmp = line.replace('"','').split(",")
        if len(tmp)==3:
            id1,desc,id2 = tmp
            if re.match('^s[35]{2}$', desc):
                if not stackings:
                    continue
                desc = "stacking_"+desc
            elif re.match('^[<>]{2}$', desc): # for use with moderna
                if not stackings:
                    continue
                desc = desc
            elif re.match('^[ct][WHS]{2}$', desc):
                if not base_pairs:
                    continue
            elif re.match('^.*(BR|BPh)$', desc):
                for new_id,old_id in [('0BPh','H_0BPh'),
                    ('1BPh','S_1BPh'),('2BPh','SW_2BPh'),
                    ('3BPh','W_3BPh'),('4BPh','W_4BPh'),
                    ('5BPh','W_5BPh'),('6BPh','W_6BPh'),
                    ('7BPh','H_7BPh'),('8BPh','H_8BPh'),
                    ('9BPh','H_9BPh'),('0BR','H_0BR'),
                    ('1BR','S_1BR'),('2BR','SW_2BR'),
                    ('3BR','W_3BR'),('4BR','W_4BR'),
                    ('5BR','W_5BR'),('6BR','W_6BR'),
                    ('7BR','H_7BR'),('8BR','H_8BR'),
                    ('9BR','H_9BR')]:
                    if new_id==desc:
                        desc = old_id
                        break
            (model1,chain1,num1,r1) = id1.split("_")[0:4]
            (model2,chain2,num2,r2) = id2.split("_")[0:4]
            if model1!='1' or model2!='1':
                continue
            if use_chain:
                name1 = chain1+num1
                name2 = chain2+num2
            else:
                name1 = num1
                name2 = num2
            res.append([(name1,name2),(r1,r2),desc])
    return res


def structure2string(structure):
    output = io.StringIO()
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output)
    return output.getvalue()

##################

def confusion_matrix_params(m):
    result = {}
    for k,v in list(m.items()):
        result[k] = v
    result['p'] = result['tp']+result['fn']
    result['n'] = result['fp']+result['tn']
    result['pp'] = result['tp']+result['fp']
    result['nn'] = result['fn']+result['tn']
    result['tpr'] = float(result['tp'])/max(1,result['tp']+result['fn'])
    result['fpr'] = float(result['fp'])/max(1,result['fp']+result['tn'])
    result['acc'] = float(result['tp']+result['tn'])/max(1,result['tp']+result['tn']+result['fp']+result['fn'])
    result['spc'] = float(result['tn'])/max(1,result['fp']+result['tn'])
    result['ppv'] = float(result['tp'])/max(1,float(result['tp']+result['fp']))
    result['npv'] = float(result['tn'])/max(1,float(result['tn']+result['fn']))
    result['fdr'] = float(result['fp'])/max(1,float(result['fp']+result['tp']))
    result['mcc'] = float(result['tp']*result['tn']-result['fp']*result['fn'])/max(1,math.sqrt(result['p']*result['n']*result['pp']*result['nn']))
    result['f1'] = float(2*result['tp'])/max(1,(result['p']+result['pp']))
    return result

##################

def compute_close_doublets(residues, min_distance=0.0, max_distance=4.0):
    result = []
    n = len(residues)
    residues_atoms = []
    residues_coords = []
    residues_num = []
    for i in range(n):
        a_dict = {}
        if residues[i] is not None:
            for a in residues[i]:
                a_dict[a.id] = a.get_coord()
        p = None
        for a in ["C1'","C2","C4","C6","P"]:
            if a in a_dict:
                p = a_dict[a]
                break
        if p is not None:
            residues_num.append(i)
            residues_atoms.append(a_dict)
            residues_coords.append(p)
    if len(residues_coords)==0:
      return []
    t = KDTree(residues_coords)
    for ii in range(len(residues_num)):
        i = residues_num[ii]
        if residues[i] is None:
            continue
        points = t.query_ball_point(residues_coords[ii], r=max_distance+15.0)
        for jj in points:
            if jj>=len(residues_num):
                continue
            j = residues_num[jj]
            if j<=i or j>=len(residues_atoms) or residues[j] is None:
                continue
            d = 1000
            for a1,p1 in list(residues_atoms[ii].items()):
                for a2,p2 in list(residues_atoms[jj].items()):
                    d = min(d, dist(p1,p2))
            if d>=min_distance and d<=max_distance:
                result.append((i,j))
    return result

##################

class DoubletDescTool(object):

    def __init__(self, description_dict):
        self.desc_dict = description_dict

    @staticmethod
    def reverse_d_id(d_id):
        tmp = d_id.split(":")
        assert len(tmp) in [2,3]
        if len(tmp)==2:
            return tmp[1]+":"+tmp[0]
        else:
            return tmp[0]+":"+tmp[2]+":"+tmp[1]

    @staticmethod
    def reverse_n_type(desc):
        return desc[::-1]

    @staticmethod
    def reverse_desc(desc):
        prefix = ""
        if len(desc)>0 and desc[0] in ['?','!']:
            prefix = desc[0]
            desc = desc[1:]
        if re.match("^[sSWH]{2}_(cis|tran)$",desc):
            return prefix+desc[1]+desc[0]+desc[2:]
        if desc=='<<':
            return prefix+'>>'
        elif desc=='>>':
            return prefix+'<<'
        elif desc in ['<>','><','']:
            return prefix+desc
        elif desc in ['diagonal-nc-ww','long-stacking-c']:
            return prefix+desc
        else:
            return prefix+'UNK_SHORT_DESC'

    @staticmethod
    def get_desc_category(desc,detailed=False,n_type=None):
        if desc is None or desc=='':
            return 'unclassified'
        if len(desc)>1 and desc[0] in ['?','!']:
            desc = desc[1:]
        if re.match("^[sSWH]{1,2}_(cis|tran)", desc):
            if detailed:
                if desc=='WW_cis' and n_type in ['CG','GC','AU','UA']:
                    return "bp-classic"
                else:
                    return "bp-non-classic"
            else:
                return "bp"
        elif re.match("^[<>]{2}", desc):
            return "stacking"
        elif re.match("^[SWH]{1,2}_\d+BP", desc):
            return "base-phosphate"
        elif re.match("^[SWH]{1,2}_\d+BR", desc):
            return "base-ribose"
        elif re.match('^(diagonal|long-stacking|long-diagonal|)', desc):
            return "other"
        elif re.match('^(base-ribose-stacking|ribose-ribose)', desc):
            return "other2"
        elif re.match('^(phosphate-stacking)', desc):
            return "other3"
        else:
            return "unknown"

    @staticmethod
    def n_type_sort_key(n_type):
        return min(n_type)+max(n_type)+n_type

    @staticmethod
    def desc_sort_key(desc):
        res = ""
        c = DoubletDescTool.get_desc_category(desc)
        if c=="bp":
            x = desc[0:2].replace("W","0").replace("H","1").replace("S","2")
            res += "1"+min(x)+max(x)+desc[3:]
        elif c=="stacking":
            res += "2"
        elif c=="base-phosphate":
            res += "3"+re.sub("^[A-Z]_","",desc)
        elif c=="base-ribose":
            res += "4"+re.sub("^[A-Z]_","",desc)
        elif c=="other" or c=="other2":
            res += "5"
        elif c=="other3":
            res += "6"
        res += desc
        return res

    def interpret_other_results(self, d):
        desc_rv = self.desc_dict['RV'].get(d.get('desc_RV',''),'')
        desc_mc = self.desc_dict['MC'].get(re.sub(r"_[0-9]+$","",d.get('desc_MC','')),'')
        desc_mo = self.desc_dict['MO'].get(d.get('desc_MO',''),'')
        desc_fr = self.desc_dict['FR'].get(d.get('desc_FR',''),'')
        if re.match('^stacking',d.get('desc_FR','')) and desc_rv=='' and desc_mc=='' and desc_mo=='':
            desc_fr=''

        if DoubletDescTool.get_desc_category(desc_rv)=='bp' and desc_rv == desc_mc:
            return ("valid", self.get_desc_category(desc_rv), desc_rv)
        elif DoubletDescTool.get_desc_category(desc_mo)=='stacking' and desc_mc == desc_mo:
            return ("valid", self.get_desc_category(desc_mo), desc_mo)
        elif DoubletDescTool.get_desc_category(desc_fr)=='base-phosphate':
            return ("valid", self.get_desc_category(desc_fr), desc_fr)
        elif DoubletDescTool.get_desc_category(desc_fr)=='base-ribose':
            return ("valid", self.get_desc_category(desc_fr), desc_fr)

        if desc_rv != '' or desc_mc != '' or desc_mo != '' or desc_fr != '':
            categories = set()
            descriptions = set()
            for d in (desc_rv,desc_mc,desc_mo,desc_fr):
                if d!='':
                    categories.add(DoubletDescTool.get_desc_category(d))
                    descriptions.add(d)
            return ('fuzzy', categories, descriptions)
        else:
            return ('unclassified','','')

def reverse_desc(desc):
    return DoubletDescTool.reverse_desc(desc)

class GraphTool(object):

    def __init__(self, filename=None, edge_type='contact',accept_fuzzy=True,filter_doublets=None):
        self.contacts = {}
        if filename is not None:
            self.load_graph(filename, edge_type, accept_fuzzy,filter_doublets)

    @staticmethod
    def _graph_to_contacts(g,accepted_type='contact',accept_fuzzy=True,filter_doublets=None):
        tmp_contacts = {}
        for (source,target,data) in g.edges(data=True):
            if data['type']==accepted_type:
                _id = source+":"+target
                desc = data['desc']
                if not accept_fuzzy and re.match(r'^\?',desc):
                    continue
                if _id not in tmp_contacts:
                    tmp_contacts[_id] = {}
                tmp_contacts[_id][desc] = data
        contacts = {}
        for _id,elems_dict in list(tmp_contacts.items()):
            contacts[_id] = list(elems_dict.values())
        return contacts

    def load_graph(self, filename, edge_type='contact', accept_fuzzy=True, filter_doublets=None):
        self.g = load_graph(filename)
        if filter_doublets is not None:
            for (source,target,data) in self.g.edges(data=True):
                _id = source+":"+target
                if _id not in filter_doublets:
                    self.g.remove_edge(source,target)

        self.contacts = GraphTool._graph_to_contacts(self.g, edge_type, accept_fuzzy, filter_doublets)

    def graph(self):
        return self.g

    def edges(self,data=False):
        return self.g.edges(data=data)

    def nodes(self,data=False):
        return self.g.nodes(data=data)

    def get_ids(self,cat=None):
        res = sorted(self.contacts.keys())
        if cat is not None:
            tmp_res = []
            for r in res:
                if self.get_contact_by_id(r, cat) != '':
                    tmp_res.append(r)
            res = tmp_res
        return res

    def get_unique_ids(self,cat=None):
        res = []
        all_ids = set(self.get_ids(cat))
        for r in all_ids:
            tmp = r.split(":")
            assert len(tmp)==2
            (chain1,num1,chain2,num2) = (tmp[0][0], int(tmp[0][1:]), tmp[1][0], int(tmp[1][1:]))
            if chain1<chain2 or (chain1==chain2 and num1<num2):
                res.append(r)
        return res

    def get_all_contacts_by_id(self,contact_id, cat=None, data=False):
        assert cat in [None,'bp','bp-classic','bp-non-classic','stacking','base-phosphate','base-ribose','other','other2','other3']
        res = self.contacts.get(contact_id, [])
        if cat is not None:
            tmp_res = []
            for r in res:
                if DoubletDescTool.get_desc_category(r['desc'],detailed=(cat in ['bp-classic','bp-non-classic']),n_type=r.get('n_type')) == cat:
                    tmp_res.append(r)
            res = tmp_res
        res.sort(key=lambda x: x['desc'])
        if data:
            return res
        else:
            return [r['desc'] for r in res]

    def get_contact_by_id(self,contact_id, cat=None, data=False):
        all_res = self.get_all_contacts_by_id(contact_id, cat, data)
        if len(all_res)==0:
            if data:
                return {'desc':'','full_desc':'','type':'contact','n_type':'??','weight':0}
            else:
                return ''
        else:
            return all_res[0]

###########################################

class GroupsTool(object):

    def __init__(self, filename=None, filter_doublets=None):
        if filename is not None:
            gr = load_json(filename)
        else:
            gr = {}
        self.fuzzy = {}
        self.ref = {}
        self.all = set()
        self.n_types = {}
        self.recognized_by = {}
        for group_name, group_elems in list(gr.items()):
            if filter_doublets is not None:
                group_elems = [x for x in group_elems if re.sub("^[^:]*:","",x) in filter_doublets]
            if re.match("^classifier/",group_name):
                (_cl, sub_category, desc, n_type) = group_name.split("/")

                rev_n_type = n_type[::-1]
                rev_desc = DoubletDescTool.reverse_desc(desc)
                rev_group_elems = []
                if sub_category in ['bp','stacking'] and n_type==rev_n_type and desc==rev_desc:
                    rev_group_elems = [GroupsTool.reverse_id(i) for i in group_elems]

                if sub_category not in self.ref:
                    self.ref[sub_category] = {}
                    if sub_category=='bp':
                        self.ref['bp-classic'] = {}
                        self.ref['bp-non-classic'] = {}

                for _id in group_elems+rev_group_elems:
                    # _id = re.sub("^([^:]):([^:]):([^:])$",r"\2:\3",full_id)
                    if _id in self.ref[sub_category] and desc!=self.ref[sub_category][_id]:
                        print("reference description conflict: %s, %s <-> %s" % (_id, desc, self.ref[sub_category][_id]))
                        continue
                    # assert self.ref[sub_category].has_key(_id)==False
                    self.ref[sub_category][_id]=desc
                    if sub_category=='bp':
                        if desc=='WW_cis' and n_type in ['CG','GC','AU','UA']:
                            self.ref['bp-classic'][_id]=desc
                        else:
                            self.ref['bp-non-classic'][_id]=desc
            if re.match("^fuzzy/",group_name):
                (_cl, sub_category, desc, n_type) = group_name.split("/")

                rev_n_type = n_type[::-1]
                rev_desc = DoubletDescTool.reverse_desc(desc)
                rev_group_elems = []
                if sub_category in ['bp','stacking'] and n_type==rev_n_type and desc==rev_desc:
                    rev_group_elems = [GroupsTool.reverse_id(i) for i in group_elems]

                if sub_category not in self.fuzzy:
                    self.fuzzy[sub_category] = {}
                for _id in group_elems+rev_group_elems:
                    # _id = re.sub("^([^:]):([^:]):([^:])$",r"\2:\3",full_id)
                    if _id not in self.fuzzy[sub_category]:
                        self.fuzzy[sub_category][_id]=[]
                    self.fuzzy[sub_category][_id].append(desc)
            if re.match("^all/all/all$",group_name):
                for _id in group_elems:
                    self.all.add(_id)
            if re.match("^(all/all|classifier/[^/]*/[^/]*)/[ACGU]{2}$",group_name):
                n_type = group_name.split("/")[-1]
                for _id in group_elems:
                    self.n_types[_id] = n_type
            if re.match("^recognized/([A-Z])+/all$",group_name):
                self.recognized_by[group_name.split("/")[1]] = set(group_elems)

    @staticmethod
    def is_reversed_id(doublet_id):
        tmp = doublet_id.split(":")
        assert len(tmp) in [2,3]
        if len(tmp)==2:
            (chain1,num1,chain2,num2) = (tmp[0][0], int(tmp[0][1:]), tmp[1][0], int(tmp[1][1:]))
        else:
            (chain1,num1,chain2,num2) = (tmp[1][0], int(tmp[1][1:]), tmp[2][0], int(tmp[2][1:]))

        if chain1<chain2 or (chain1==chain2 and num1<num2):
            return False
        else:
            return True

    @staticmethod
    def reverse_id(_id):
        tmp = _id.split(":")
        if len(tmp)==2:
            return tmp[1]+":"+tmp[0]
        else:
            return tmp[0]+":"+tmp[2]+":"+tmp[1]

    def get_all_ids(self):
        return self.all

    def get_n_type(self,doublet_id):
        if doublet_id in self.n_types:
            return self.n_types[doublet_id]
        rev_doublet_id = re.sub("^([^:]+):([^:]+):([^:]+)$",r"\1:\3:\2", doublet_id)
        if rev_doublet_id in self.n_types:
            return self.n_types[rev_doublet_id][::-1]
        return '??'

    def get_ref_ids(self,sub_category):
        if sub_category not in self.ref:
            return []
        return list(self.ref[sub_category].keys())

    def get_fuzzy_ids(self,sub_category):
        if sub_category not in self.fuzzy:
            return []
        return list(self.fuzzy[sub_category].keys())

    def get_ref_desc(self,doublet_id,sub_category):
        if sub_category in ['bp-classic','bp-non-classic']:
            tmp_res = self.get_ref_desc(doublet_id,'bp')
            n_type = self.get_n_type(doublet_id)
            is_classic = (tmp_res=='WW_cis' and n_type in ['CG','GC','AU','UA'])
            if sub_category=='bp-classic' and is_classic:
                return tmp_res
            elif sub_category=='bp-non-classic' and not is_classic:
                return tmp_res
            else:
                return ''
        if sub_category not in self.ref:
            return ''
        return self.ref[sub_category].get(doublet_id,'')

    def get_fuzzy_desc(self,doublet_id,sub_category):
        if sub_category not in self.fuzzy:
            return []
        return self.fuzzy[sub_category].get(doublet_id,[])

    def is_valid(self,doublet_id):
        return doublet_id in self.all

    def get_recognized_by(self,prg):
        if prg not in self.recognized_by:
            return set([])
        else:
            return self.recognized_by[prg]

#########################

class ResiduesDict(object):

    def __init__(self,data_dir=DATA_DIR,reduced_atoms=None):
        self.full_atoms = False
        if reduced_atoms is not None:
            if '*' in reduced_atoms:
                self.full_atoms = True
            self.reduced_atoms = reduced_atoms
        else:
            self.reduced_atoms = ('C2','C4','C6',"C1'","O2'","C2'","P","OP1","OP2")
        self.data_dir = data_dir
        self.residues = {}
        self.residues_name = {}
        self.residues_org_id = {}

    def load_pdb(self,pdbid,residues=None):
        fn = PDBObject.pdb_fn(pdbid,"residues",data_dir=self.data_dir)
        json = load_json(fn)
        if self.full_atoms:
            data = dict([(k,v['atoms']) for k,v in list(json.items())])
        else:
            data = dict([(k,self._reduce_atoms(v['atoms'])) for k,v in list(json.items())])
        data_name = dict([(k,v['resname']) for k,v in list(json.items())])
        data_org_id = dict([(k,v.get('org_id')) for k,v in list(json.items())])
        json = None
        if residues is not None:
            for k in list(data.keys()):
                if k not in residues:
                    del data[k]
                    del data_name[k]
                    del data_org_id[k]
        self.residues[pdbid] = data
        self.residues_name[pdbid] = data_name
        self.residues_org_id[pdbid] = data_org_id

    def _reduce_atoms(self,atoms):
        return dict([(k,v) for k,v in list(atoms.items()) if k in self.reduced_atoms])

    def get(self,rid):
        (pdbid,num) = rid.split(":")
        if pdbid not in self.residues:
            self.load_pdb(pdbid)
        res = self.residues[pdbid].get(num)
        return res

    def get_resname(self,rid):
        (pdbid,num) = rid.split(":")
        if pdbid not in self.residues_name:
            self.load_pdb(pdbid)
        return self.residues_name[pdbid].get(num)

    def get_org_id(self,rid,full_id=True):
        (pdbid,num) = rid.split(":")
        if pdbid not in self.residues_org_id:
            self.load_pdb(pdbid)
        r = self.residues_org_id[pdbid].get(num,None)
        if r is not None:
            if full_id:
                return pdbid+":"+r
            else:
                return r
        else:
            return None

class DoubletsDict(object):

    def __init__(self,data_dir=DATA_DIR,reduced_atoms=None):
        self.data_dir = data_dir
        self.rd = ResiduesDict(data_dir,reduced_atoms)

    def __str__(self):
        return "DoubletsDict()"

    def load_pdb(self,pdbid,doublets=None):
        residues = None
        if doublets is not None:
            residues = set([x.split(":")[1] for x in doublets]+[x.split(":")[2] for x in doublets])
        self.rd.load_pdb(pdbid,residues)

    def load_pdb_files(self,ids,verbose=False):
        from .progressbar import ProgressBar, Percentage, Bar, ETA

        all_pdbs = {}
        for d_id in ids:
            pdb_id = d_id.split(":")[0]
            if pdb_id not in all_pdbs:
                all_pdbs[pdb_id]=[]
            all_pdbs[pdb_id].append(d_id)
        if len(all_pdbs)==0:
            return
        if verbose:
            widgets = ['load residues',' ', Percentage(), ' ', Bar(), ' ', ETA()]
            pbar = ProgressBar(widgets=widgets, maxval=len(list(all_pdbs.keys()))).start()
        for i,(pdb_id,doublets) in enumerate(all_pdbs.items()):
            self.load_pdb(pdb_id,doublets)
            if verbose:
                pbar.update(i)
        if verbose:
            pbar.finish()

    def get(self,did):
        (pdbid,num1,num2) = did.split(":")
        r1 = self.rd.get(pdbid+":"+num1)
        r2 = self.rd.get(pdbid+":"+num2)
        return (r1,r2)

    def get_normalized(self,did):
        n_type = self.get_n_type(did)
        (r1,r2) = self.get(did)
        return normalize_points((r1,r2),n_type[0].upper())

    def get_fit(self,did):
        n_type = self.get_n_type(did)
        (r1,r2) = self.get(did)
        return fit_points((r1,r2),n_type)

    def get_n_type(self,did):
        (pdbid,num1,num2) = did.split(":")
        r1 = self.rd.get_resname(pdbid+":"+num1)
        r2 = self.rd.get_resname(pdbid+":"+num2)
        if r1 is None or r2 is None:
            return None
        return r1+r2

    def get_residue(self,did):
        return self.rd.get(did)

    def get_residue_name(self,did):
        return self.rd.get_resname(did)

    def get_org_id(self,did):
        (pdbid,num1,num2) = did.split(":")
        r1 = self.rd.get_org_id(pdbid+":"+num1,full_id=False)
        r2 = self.rd.get_org_id(pdbid+":"+num2,full_id=False)
        if r1 is None or r2 is None:
            return None
        return pdbid+":"+r1+":"+r2

    def get_org_ids(self,did_iter):
        return [self.get_org_id(did) for did in did_iter]

#################################

class FileNamesObject:

    @staticmethod
    def groups_fn(setname=None,reduced=None,cross_validation_num=None,pdb_id=None,data_dir=DATA_DIR):
        assert setname in [None,'all','bench','training']

        if pdb_id is not None:
            return PDBObject.pdb_fn(pdb_id,"groups",data_dir)

        if setname is not None and reduced is not None:
            if reduced:
                p1 = "reduced_groups"
            else:
                p1 = "groups"
            if cross_validation_num is not None:
                p2 = "-cross%s" % cross_validation_num
                if setname=="training":
                    p2 += "t"
                else:
                    p2 += "v"
            else:
                if setname=='training':
                    p2 = ""
                else:
                    p2 = "-%s" % setname
            return os.path.join(data_dir,"groups","rna."+p1+p2+".json.gz")
        raise Exception("Unknown file type: %s" % locals())

    @staticmethod
    def eval_fn(setname=None,t=None,pdb_id=None,data_dir=DATA_DIR):
        if pdb_id is not None:
            assert t in ['eval1','eval2','eval3','corr']
            return PDBObject.pdb_fn(pdb_id,t,data_dir)

        assert setname in ['small','medium','training','bench','all']

        if t=='eval1':
            return os.path.join(data_dir,"eval","new-classifier-eval-"+setname+".json.gz")
        if t=='eval2':
            return os.path.join(data_dir,"eval","new-classifier-full-stats-"+setname+".json.gz")
        if t=='eval3':
            return os.path.join(data_dir,"eval","new-classifier-full-stats2-"+setname+".json.gz")
        if t=='corr':
            return os.path.join(data_dir,"eval","new-classifier-corr-"+setname+".json.gz")
        if t in ['eval1_in','eval2_in','eval3_in','corr_in']:
            return os.path.join(data_dir,"tmp","eval_in","new-classifier-corr-"+setname+"-"+t+".in")
        raise Exception("Unknown file type: %s" % locals())


class PDBObject:

    def __init__(self,pdb_id,data_dir=DATA_DIR):
        self.data_dir = data_dir
        self.pdb_id = pdb_id.lower()

    @staticmethod
    def pdb_fn(pdb_id,t,data_dir=DATA_DIR):
        pdb_id = pdb_id.lower()
        if t=='pdb':
            return os.path.join(data_dir,"pdb_files",pdb_id+"_rna.pdb.gz")
        elif t=='residues':
            return os.path.join(data_dir,"residues",pdb_id+"_rna.residues.json.gz")
        elif t=='rev_desc':
            return os.path.join(data_dir,"tmp","rev_desc",pdb_id+"_rna.rev_desc.json.gz")
        elif t=='res_map':
            return os.path.join(data_dir,"pdb_files",pdb_id+"_rna.json.gz")
        elif t=='groups':
            return os.path.join(data_dir,"groups_by_pdb",pdb_id+"_rna.groups.json.gz")
        elif t=='contacts_CL':
            return os.path.join(data_dir,"contacts",pdb_id+"_rna.new-classifier-graph.json.gz")
        elif t=='contacts_FR':
            return os.path.join(data_dir,"contacts",pdb_id+"_rna.contacts.fr3d-graph.json.gz")
        elif t=='contacts_MC':
            return os.path.join(data_dir,"contacts",pdb_id+"_rna.contacts.mc-annotate-graph.json.gz")
        elif t=='contacts_MO':
            return os.path.join(data_dir,"contacts",pdb_id+"_rna.contacts.moderna-graph.json.gz")
        elif t=='contacts_RV':
            return os.path.join(data_dir,"contacts",pdb_id+"_rna.contacts.rnaview-graph.json.gz")
        elif t=='close_doublets':
            return os.path.join(data_dir,"contacts",pdb_id+"_rna.close-doublets-graph.json.gz")
        elif t=='eval1':
            return os.path.join(data_dir,"tmp","eval",pdb_id+"_rna.new-classifier-eval.json.gz")
        elif t=='eval2':
            return os.path.join(data_dir,"tmp","eval",pdb_id+"_rna.new-classifier-stats.json.gz")
        elif t=='eval3':
            return os.path.join(data_dir,"tmp","eval",pdb_id+"_rna.new-classifier-stats2.json.gz")
        elif t=='corr':
            return os.path.join(data_dir,"tmp","corr",pdb_id+"_rna.classifiers-corr.json.gz")
        else:
            raise Exception("Unknown file type: %s" % t)

    def fn(self,t):
        return PDBObject.pdb_fn(self.pdb_id,t,data_dir=self.data_dir)
