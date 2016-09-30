#!/usr/bin/env python
#
# Caution!
# this script requires Residues Dictionaries in gc-data/residues/
# (for example: gc-data/residues/1rna_rna.residues.json.gz)
#
# Example usage:
#
# - generate single doublet
#   ./gen-pdb.py --doublet-id=157D:A10:B1 -o example-ww-cis.pdb --type=full
#
# - generate multiple doublets from cmd line
#   ./gen-pdb.py --doublet-ids=157D:A10:B1,157D:A2:B9,157D:B10:A1 -o example-ww-cis2.pdb --type=bp
#
# - generate multiple doublets from JSON list
#   echo '["157D:A10:B1", "157D:A2:B9", "157D:B10:A1", "157D:B2:A9"]' > a.json
#   ./gen-pdb.py --input-json=a.json -o example-ww-cis3.pdb --type=full 
#
# - generate multiple doublets from JSON dict
#   ./gen-pdb.py --input-json=gc-data/groups/rna.reduced_groups.json.gz --filter-keys=classifier/bp/WW_cis/CG \
#     -o example-ww-cis3.pdb.gz --type=full --limit=500
#   ./gen-pdb.py --input-json=gc-data/groups/rna.reduced_groups.json.gz --filter-keys='classifier/stacking/>>/AU' \
#     -o example-stacking.pdb.gz --type=full --limit=500
import re
import sys
import itertools
import multiprocessing
import random
from optparse import OptionParser

from Bio.SVDSuperimposer import SVDSuperimposer

import scipy as scipy
import numpy as np
import scipy.cluster.hierarchy as sch
from utils import *
from distances import rmsd_distance, doublet_params_dict, normalize_points, center_vector, \
    ribose_center_vector, NORMALIZED_BASE
from progressbar import ProgressBar, Percentage, Bar, ETA

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""tool for generating PDB files with doublets""")
    parser.add_option("--doublet-id", dest="doublet_id",
                  help="doublet id", metavar="PDBID:R1:R2")
    parser.add_option("--doublet-ids", dest="doublet_ids",
                  help="doublet id list (comma separated)", metavar="PDBIDS")
    parser.add_option("--input-json", dest="input_json",
                  help="read doublets from json file", metavar="FILE")
    parser.add_option("--group-info", dest="group_info",
                  help="read doublets from group info file", metavar="FILE")
    parser.add_option("--filter-keys", dest="filter_keys",
                  help="comma separated list of keys", metavar="KEYS")
    parser.add_option("--type", dest="type",
                  help="base type (bp,bp-centers,stacking,base-ribose,base-phosphate,full)", metavar="T", default="bp")
    parser.add_option("--only-chain-b", dest="only_chain_b", action="store_true",
                  help="save only chain B", default=False)
    parser.add_option("--add-normalized-base", dest="add_normalized_base", action="store_true",
                  help="add normalized base", default=False)
    parser.add_option("--add-stacking-convex", dest="add_stacking_convex", action="store_true",
                  help="add stacking convex", default=False)
    parser.add_option("-o","--output", dest="output",
                  help="save output to file", metavar="FILE")
    parser.add_option("--data-dir", dest="data_dir",
                  help="data directory", metavar="DIR", default="gc-data")
    parser.add_option("--limit", dest="limit",
                  help="models limit", metavar="N")
    parser.add_option("--randomize", dest="randomize", action='store_true',
                  help="randomize models")
    parser.add_option("--save-json", dest="save_json", action='store_true',
                  help="save json file")
    parser.add_option("--translate-doublet-ids", dest="translate_doublet_ids", action='store_true',
                  help="save json file translated doublet ids (to original names)")
    parser.add_option("--dont-normalize", dest="dont_normalize", action="store_true",
                  help="dont normalize location", default=False)
    parser.add_option("--skip-empty", dest="skip_empty", action="store_true",
                  help="skip generation of empty files", default=False)
    parser.add_option("--quiet", dest="quiet", action="store_true",
                  help="quiet version", default=False)
    (options, args)  = parser.parse_args()
    return (parser, options, args)

def points2model(model_num,points,n_type,first_atom_num=0):
    num = first_atom_num
    model = PDB.Model.Model(model_num)
    for i,p in enumerate(points):
        chain = PDB.Chain.Chain(chr(ord('A')+i))
        num += 1
        residue = PDB.Residue.Residue((' ',i,' '),n_type[i],num)
        for j,(k,v) in enumerate(p.items()):
            num += 1
            kk = k
            if "NEXT:" in kk:
                kk = k.replace("NEXT:","")
                if p.has_key(kk):
                    kk += "N"
            element = k.replace("NEXT:","").replace("1","").replace("2","")[0]
            # print "serial: %s, k=%s, element=%s" % (num,k,element)
            atom = PDB.Atom.Atom(kk,v,1,1," ",fullname=kk,serial_number=num,element=element)
            residue.add(atom)
        chain.add(residue)
        model.add(chain)
    return model
    
def stacking_convex_points(n_type):
    points = {
        'A': zip(
                [ -2.100463,  0.000000,  4.271447,  1.920945,  0.230436, -2.100463],
                [  0.447145, -1.009320,  1.317924,  5.150733,  4.699718,  0.447145]
             ),
        'C': zip(
                [ -2.082733,  0.000000,  2.269450,  1.203833, -0.527970, -2.036772, -2.082733],
                [  0.123632, -1.010259, -0.120783,  4.411996,  4.602202,  2.647095,  0.123632]
             ),
        'G': zip(
                [ -2.101572,  0.000000,  4.872516,  5.295175,  3.613335,  1.396986, -0.751391, -2.101572],
                [  0.463584, -1.009529,  0.097781,  1.782283,  3.374547,  4.394213,  2.120132,  0.463584]
             ),
        'U': zip(
                [ -2.082780,  0.000000,  2.292490,  2.092152,  0.177156, -2.124577, -2.082780],
                [  0.111836, -1.008947, -0.048394,  2.445179,  4.020060,  2.616537,  0.111836]
             ),
    }
    return dict( [ ('P%d'%i,(x,y,0.0)) for i,(x,y) in enumerate(points[n_type[0]]) ] )


def gen_pdb(doublets,options):
    limit = None
    if options.limit:
        limit = int(options.limit)

    base_atoms = set(sum([NORMALIZED_BASE[c].keys() for c in ['A','C','G','U']],[]))
    ribose_atoms = set(["O2'","O3'","O4'","C1'","C2'","C3'","C4'",])
    phosphate_atoms = set(["OP1","OP2","O5'","NEXT:O3'","P"])

    r_atoms = list(base_atoms)
    filter0 = base_atoms
    filter1 = base_atoms
    if options.type in ['base-ribose']:
        r_atoms += list(ribose_atoms)
        filter1 = ribose_atoms
    elif options.type in ['base-phosphate']:
        r_atoms += list(phosphate_atoms)
        filter1 = phosphate_atoms
    if options.type in ['ribose-ribose']:
        r_atoms += list(ribose_atoms)
        filter0 = ribose_atoms
        filter1 = ribose_atoms
    elif options.type in ['full']:
        filter0 = None
        filter1 = None
        r_atoms = ['*']
    elif options.type=='bp-centers':
        filter0 = base_atoms
        filter1 = ['CC']
    elif options.type=='ribose-centers':
        r_atoms += list(ribose_atoms)
        filter0 = base_atoms
        filter1 = ['RC']

    if limit is not None:
        doublets = doublets[0:limit]

    dd = DoubletsDict(options.data_dir, reduced_atoms=r_atoms)
    dd.load_pdb_files(doublets, verbose=not options.quiet)
    
    model_num = 0
    atom_num = 0
    s = PDB.Structure.Structure("superimpose")
    if options.add_normalized_base and len(doublets)>0:
        n_type = dd.get_n_type(doublets[0])
        n_type0 = n_type[0]
        model_num += 1
        model = points2model(model_num,(NORMALIZED_BASE[n_type0],{}),n_type,first_atom_num=atom_num)
        s.add(model)
        atom_num += len([x for x in model.get_atoms()])
    if options.add_stacking_convex and len(doublets)>0:
        n_type = dd.get_n_type(doublets[0])
        model_num += 1
        model = points2model(model_num,(stacking_convex_points(n_type),{}),n_type,first_atom_num=atom_num)
        s.add(model)
        atom_num += len([x for x in model.get_atoms()])
    saved_doublets = []
    for d_id in doublets:
        n_type = dd.get_n_type(d_id)
        if n_type is None:
            print "WARNING! Unknown N_TYPE for %s, ignoring doublet" % d_id
            continue
        if options.dont_normalize:
            new_p = dd.get(d_id)
        else:
            new_p = dd.get_normalized(d_id)
        # add extra points:
        for i,fset in (0,filter0),(1,filter1):
            if fset is not None and 'CC' in fset:
                new_p[i]['CC'] = center_vector(new_p[i],n_type[1])
            if fset is not None and 'RC' in fset:
                new_p[i]['RC'] = ribose_center_vector(new_p[i])
        if filter0 is not None:
            for k,v in new_p[0].items():
                if k not in filter0:
                    del new_p[0][k]
        if filter1 is not None:
            for k,v in new_p[1].items():
                if k not in filter1:
                    del new_p[1][k]
        model_num += 1
        if limit is not None and model_num>limit:
            break
        
        if not options.quiet:
            print "adding %s (%s) as model %d" % (d_id, n_type, model_num)
        if options.only_chain_b:
            new_p[0] = {}

        model = points2model(model_num,new_p,n_type,first_atom_num=atom_num)
        s.add(model)
        atom_num += len([x for x in model.get_atoms()])
        saved_doublets.append(d_id)
    out_fn = options.output
    save_pdb(out_fn,s)
    if options.save_json:
        fn = re.sub(".pdb(.gz)?$","",out_fn) + ".json"
        if options.translate_doublet_ids:
            saved_doublets = [dd.get_org_id(x) for x in saved_doublets]
        save_json(fn,saved_doublets)
        

########################

def main():
    (parser,options,_args) = parse_args()
    
    doublets = []
    if options.doublet_id:
        doublets = [options.doublet_id]
    elif options.doublet_ids:
        doublets = options.doublet_ids.replace('"',"").replace(" ","").split(",")
    elif options.input_json:
        json = load_json(options.input_json)
        if isinstance(json,list):
            doublets = json
        elif isinstance(json,dict):
            if options.filter_keys is None:
                filter_keys = None
            else:
                filter_keys = options.filter_keys.split(",")
            if filter_keys is not None and len(filter_keys)==1:
                regexp = re.compile("^"+options.filter_keys+"$")
                for key,values in json.items():
                    if not regexp.match(key):
                        continue
                    doublets += values
            else:
                for key,values in json.items():
                    if (filter_keys is not None) and (not (key in filter_keys)):
                        continue
                    doublets += values
            doublets = sorted(list(set(doublets)))
        else:
            raise Exception("Unknown json format")
    elif options.group_info:
        json = load_json(options.group_info)
        assert isinstance(json,dict)
        assert len(json.keys())==1
        key = json.keys()[0]
        v = json.values()[0]
        options.type = key.split("/")[1]
        if options.filter_keys is None:
            options.filter_keys = "doublets"
        filter_keys = options.filter_keys.split(",")
        for key,values in v.items():
            if not (key in filter_keys):
                continue
            if key in ['doublets','all_doublets']:
                doublets += values
            elif key in ['neigh_unclassified','neigh_other']:
                doublets += [d_id for row in values for d_id,dist in row]
            else:
                raise Exception("Unknown key: %s" % key)
        doublets = sorted(list(set(doublets)))
    else:
        parser.error("Specify doublets to generate")

    if options.skip_empty and len(doublets)==0:
        print "skipping generation of empty files!"
        return
    if options.randomize:
        random.shuffle(doublets)

    gen_pdb(doublets,options)
    
if __name__=="__main__":
    main()