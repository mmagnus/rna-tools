#!/usr/bin/env python
import sys
import os
import re
import shutil
import StringIO
import simplejson
import gzip
import tempfile
import networkx as nx
import itertools
from networkx.readwrite import json_graph
from optparse import OptionParser
from normalize_pdb import normalize_pdb
from itertools import combinations
from structure_ciach import StructureCiachCiach

from utils import *
from cl_settings import *

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""extract contacts from PDB file
""")
    parser.add_option("-o", "--output-dir", dest="output_dir",
                  help="save motifs to given directory", metavar="DIR",
                  default="contacts-db")
    parser.add_option("--output-graph", dest="output_graph",
                  help="save doublets as NetworkX graph (in JSON format)", metavar="FILE")
    parser.add_option("--output-json", dest="output_json",
                  help="save motifs to given JSON file", metavar="FILE")
    parser.add_option("-i", "--input", dest="input",
                  help="parse structures from directory or use given file", metavar="DIR|FILE",
                  default="pdb-simrna")
    parser.add_option("--use-rnaview", dest="use_rnaview", action="store_true",
                  help="use rna view to detect contacts",
                  default=False)
    parser.add_option("--use-mc-annotate", dest="use_mc_annotate", action="store_true",
                  help="use MC-Annotate to detect contacts",
                  default=False)
    parser.add_option("--use-moderna", dest="use_moderna", action="store_true",
                  help="use Moderna algorithm to detect contacts",
                  default=False)
    parser.add_option("--use-new-alg2", dest="use_new_alg2", action="store_true",
                  help="use new alg2 to detect contacts",
                  default=False)
    parser.add_option("--use-new-alg2-mc", dest="use_new_alg2_mc", action="store_true",
                  help="use new alg2 to detect contacts (with mc database)",
                  default=False)
    parser.add_option("--use-new-alg3", dest="use_new_alg3", action="store_true",
                  help="use new alg3 to detect contacts (with rnaview database)",
                  default=False)
    parser.add_option("--use-fr3d", dest="use_fred", action="store_true",
                  help="use FR3D algorithm (local)",
                  default=False)
    parser.add_option("--use-rna-hub", dest="use_rna_hub", action="store_true",
                  help="use FR3D algorithm (retrive results from Web)",
                  default=False)
    parser.add_option("--negative", dest="negative", action="store_true",
                  help="also only negative matches (close pairs without contact)",
                  default=False)
    parser.add_option("--reverse-pdb", dest="reverse", action="store_true",
                  help="also generate contacts in reversed PDB files",
                  default=False)
    parser.add_option("--with-backbone", dest="with_backbone", action="store_true",
                  help="add backbone fragments of the contacts",
                  default=False)
    parser.add_option("--no-pdb", dest="no_pdb", action="store_true",
                  help="do not save PBD info to JSON files",
                  default=False)
    parser.add_option("--extract-all", dest="extract_all", action="store_true",
                  help="extract all kind of contacts",
                  default=False)
    parser.add_option("--extract-stackings", dest="extract_stackings", action="store_true",
                  help="extract stackings only",
                  default=False)
    parser.add_option("--close-doublets", dest="close_doublets", action="store_true",
                  help="extract close doublets",default=False)
    parser.add_option("--min-distance", dest="min_distance", metavar="DIST",
                  help="minimal distance",default="0.0")
    parser.add_option("--max-distance", dest="max_distance", metavar="DIST",
                  help="maximal distance",default="4.0")
    parser.add_option("--dont-normalize", dest="dont_normalize", action="store_true",
                  help="do not normalize input files",default=False)
    parser.add_option("--residue-mapping", dest="residue_mapping", metavar="FILE",
                  action="store",help="use residue mapping from FILE")
    parser.add_option("--descriptions-dict", dest="descriptions_dict", metavar="JSON_FILE",
                  help="load descriptions dictionary from JSON file")
    parser.add_option("--skip-invalid-doublets", dest="skip_invalid", action="store_true",
                  help="skip invalid doublets (with missing atoms/or too distant)",default=False)
    (options, args)  = parser.parse_args()

    if options.use_mc_annotate and options.reverse:
        options.reverse_id = True
    else:
        options.reverse_id = False
    if not options.descriptions_dict:
        options.descriptions_dict = os.path.join(CLARNA_DIR,"descriptions-dict.json")

    return (parser, options, args)

def find_close_contacts(pdb, min_distance, max_distance, options):
    result = []
    structure = load_pdb(pdb)

    if re.match('^.*/zz.*.pdb(.gz)$', options.input):
        req_atoms_list = ('C2','C4','C6')
    else:
        req_atoms_list = REQ_ATOMS_LIST

    ciach = StructureCiachCiach(structure,options.dont_normalize,req_atoms_list=req_atoms_list)
    residues = [r for r in structure[0].get_residues()]
    n = len(residues)
    r_id = []
    r_name = []
    for i,r in enumerate(residues):
        _id = ""
        if options.dont_normalize:
            _id += r.get_parent().get_id()
        _id += (str(r.get_id()[1]) + r.get_id()[2]).strip()
        _name = r.resname.strip()
        r_id.append(_id)
        r_name.append(_name)

    for (i,j) in compute_close_doublets(residues,min_distance,max_distance):
        id1,id2 = r_id[i],r_id[j]
        resname1,resname2 = r_name[i],r_name[j]
        if ciach.good_residues_num.has_key(id1) and ciach.good_residues_num.has_key(id2):
            result.append([(id1, id2), (resname1, resname2), "close-doublet"])
        else:
            print "skipping %s:%s (bad residues)" % (id1,id2)
    return result

def save_output_graph(contacts,residues,residues_info,prg,options):
    descriptions_dict = load_json(options.descriptions_dict)
    g = nx.MultiDiGraph()
    for r in residues:
        r_id = ""
        if options.dont_normalize:
            r_id += r.get_parent().get_id()
        r_id += str(r.get_id()[1])
        resname = r.resname.strip()
        kwargs = {'resname':resname}
        if residues_info.has_key(r_id):
            kwargs['conf'] = residues_info[r_id]['conf']
        g.add_node(r_id,**kwargs)
    if prg=='--':
        edge_type='dist'
    else:
        edge_type='contact'
    all_nodes = set(g.nodes())
    for (num1,num2),(r1,r2),desc in contacts:
        tmp_desc = desc
        if prg=='--':
            short_desc='close-doublet'
        else:
            if prg=='MC':
                tmp_desc = re.sub("_\d+$","",tmp_desc)
            elif prg=='FR':
                if re.match('^n',tmp_desc):
                    tmp_desc = ""
                    desc = ""
            short_desc = descriptions_dict[prg].get(tmp_desc,'UNK_SHORT_DESC')
        if prg=="RV":
            # RNA-view uppercase all chains, so we should check the lowercase version
            _num1 = num1[0].lower()+num1[1:]
            _num2 = num2[0].lower()+num2[1:]
            if num1 not in all_nodes and _num1 in all_nodes:
                num1 = _num1
            if num2 not in all_nodes and _num2 in all_nodes:
                num2 = _num2
        n_type=r1+r2
        if desc!="":
            g.add_edge(num1,num2,type=edge_type,prg=prg,desc=short_desc,full_desc=desc,n_type=n_type)
        if prg!='FR':
            if prg=='--':
                rev_short_desc = 'close-doublet'
            else:
                rev_short_desc = reverse_desc(short_desc)
            g.add_edge(num2,num1,type=edge_type,prg=prg,desc=rev_short_desc,full_desc="REV:"+desc,reverse=True,n_type=n_type[::-1])
    if re.match(r"^.*\.gz$", options.output_graph):
        f = gzip.open(options.output_graph,"w")
    else:
        f = open(options.output_graph,"w")
    json_graph.dump(g, f, indent=2)
    f.close()

def extract_contacts(pdb, id, numbers, options):
    extract_bp = not options.extract_stackings
    extract_stackings = options.extract_stackings
    if options.extract_all:
        extract_bp = True
        extract_stackings = True
        
    residues_info = {}

    if options.use_rnaview:
        prg = "RV"
        contacts = parse_rnaview(run_rnaview(pdb),
            base_pairs=extract_bp,
            stackings=extract_stackings,
            use_chain=options.dont_normalize
        )
    elif options.use_mc_annotate:
        prg = "MC"
        contacts, residues_info = parse_mc_annotate(run_mc_annotate(pdb),
            base_pairs=extract_bp,
            stackings=extract_stackings,
            use_chain=options.dont_normalize
        )
    elif options.use_moderna:
        prg = "MO"
        assert  extract_stackings==True
        res_str = run_moderna(pdb)
        contacts = parse_fr3d(res_str,
            base_pairs=extract_bp,
            stackings=extract_stackings,
            use_chain=options.dont_normalize
        )
    elif options.use_new_alg2:
        prg = "ZZ"
        contacts = parse_rnaview(run_new_alg2(pdb),my_output=True)
    elif options.use_new_alg2_mc:
        prg = "ZZ"
        contacts = parse_rnaview(run_new_alg2(pdb,db="mc-annotate"),my_output=True)
    elif options.use_new_alg3:
        prg = "ZZ"
        contacts = parse_rnaview(run_new_alg3(pdb,db="rna-view.json"),my_output=True)
    elif options.use_fred:
        prg = "FR"
        contacts = parse_fr3d(run_fr3d(pdb),
            base_pairs=extract_bp,
            stackings=extract_stackings)
    elif options.use_rna_hub:
        prg = "FR"
        contacts = parse_fr3d(download_rna_hub(id[0:4]),
            base_pairs=extract_bp,
            stackings=extract_stackings)
        # TODO use mapping!

        if options.residue_mapping:
            num_dict = {}
            for chain_id,chain_elements in load_json(options.residue_mapping).items():
                for num,r in chain_elements.items():
                    num_dict["%s%s%s" % (chain_id,r['resseq'],r['icode'])] = "%s%s" % (chain_id,num)
            new_contacts = []
            for ((num1,num2),(r1,r2),desc) in contacts:
                if not num_dict.has_key(num1):
                    print "missing info about residue: %s" % num1
                    continue
                if not num_dict.has_key(num2):
                    print "missing info about residue: %s" % num2
                    continue
                new_contacts.append(((num_dict[num1],num_dict[num2]),(r1,r2),desc))
            contacts = new_contacts
        
    elif options.close_doublets:
        prg = "--"
        min_distance = float(options.min_distance)
        max_distance = float(options.max_distance)
        contacts = find_close_contacts(pdb, min_distance, max_distance, options)
    else:
        contacts = []

    if options.skip_invalid:
        cc = find_close_contacts(pdb,0.0,4.0,options)
        valid_doublets = set([doublet_id for (doublet_id,_a,_b) in cc])
        valid_doublets = valid_doublets.union(set([doublet_id[::-1] for (doublet_id,_a,_b) in cc]))
        if prg=='RV':
            valid_doublets = valid_doublets.union(set([(d1.upper(),d2.upper()) for d1,d2 in valid_doublets]))
        org_contacts = contacts
        contacts = list(itertools.ifilter(lambda x: x[0] in valid_doublets,org_contacts))

    structure = load_pdb(pdb)
    residues = [r for r in structure[0].get_residues()]
    n = len(residues)

    if options.output_graph:
        save_output_graph(contacts,residues,residues_info,prg,options)
        return
        
    ciach = StructureCiachCiach(structure,options.dont_normalize)


    pairs = set()
    json_result = []
    for c in contacts:
        ((num1,num2),(r1,r2),desc) = c
        if r1 == "-":
            r1 = ciach.get_resname(num1)
        if r2 == "-":
            r2 = ciach.get_resname(num2)
        if not options.negative:
            if options.reverse_id:
                (pdb_num1,pdb_num2) = (num1,num2)
                (fn_num1,fn_num2) = (n-int(num1),n-int(num2))
            else:
                (fn_num1,fn_num2) = (num1,num2)
                (pdb_num1,pdb_num2) = (num1,num2)
            n_type = r1+r2
            part_fn = "%s/%s%s_desc_%s_%s_from_%s_num_%s_%s.pdb" % (options.output_dir,r1,r2,prg,desc,id,fn_num1,fn_num2)
            if options.output_json:
                if options.no_pdb:
                    mode = "dict-no-pdb"
                else:
                    mode="dict"
                d = ciach.extract(mode, pdb_num1, pdb_num2, desc, n_type, prg, options.with_backbone)
                if d:
                    print "saving %s to json" % os.path.basename(part_fn)
                    json_result.append(d)
            else:
                print "extracting %s " % part_fn
                if not ciach.extract(part_fn, pdb_num1, pdb_num2, desc, n_type, prg, options.with_backbone):
                    print "FAIL"
        pairs.add((num1,num2))
        pairs.add((num2,num1))

    if options.negative:
        residues = [(simplify_residue(r),r.get_resname().strip(), str(r.get_id()[1])) for r in structure.get_residues()]
        for (i,j) in combinations(xrange(len(residues)), 2):
            (s1,r1,num1) = residues[i]
            (s2,r2,num2) = residues[j]
            if res_c1p_distance(s1, s2) > 15:
                continue
            if (num1,num2) not in pairs:
                part_fn = "%s/%s%s_negative_%s_from_%s_num_%s_%s_.pdb" % (options.output_dir,r1,r2,prg,id,num1,num2)
                print "extracting %s " % part_fn
                if not ciach.extract(part_fn, num1, num2, "negative", n_type, '--', options.with_backbone):
                    print "FAIL"
    if options.output_json:
        if re.match(r"^.*.gz$",options.output_json):
            f = gzip.open(options.output_json,"wb")
        else:
            f = open(options.output_json,"w")
        simplejson.dump(json_result, f)

def handle_file(pdb, options):
    print "Processing file: %s" % pdb
    (_fd,tmp_pdb) = tempfile.mkstemp(".norm.pdb")
    id = os.path.basename(pdb).replace(".pdb","")
    if options.dont_normalize:
        numbers = {}
        if re.match(r"^.*.gz$", pdb):
            f1 = gzip.open(pdb,"rb")
            f2 = open(tmp_pdb,"wb")
            f2.write(f1.read())
            f1.close()
            f2.close()
        else:
            shutil.copyfile(pdb, tmp_pdb)
    else:
        numbers = normalize_pdb(pdb, tmp_pdb, reverse=options.reverse, reverse_id=options.reverse_id)
    extract_contacts(tmp_pdb, id, numbers, options)
    os.unlink(tmp_pdb)


def main():
    (parser, options, args) = parse_args()
    if (not options.use_mc_annotate and
        not options.use_rnaview and
        not options.use_moderna and
        not options.use_new_alg2 and
        not options.use_new_alg3 and
        not options.use_new_alg2_mc and
        not options.use_fred and
        not options.use_rna_hub and
        not options.close_doublets):
        print "select method"
        parser.print_help()
        exit(1)

    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)

    if os.path.isdir(options.input):
        files = pdb_files(options.input)
    elif os.path.isfile(options.input):
        files = [options.input]
    else:
        print "specify --input"
        parser.print_help()
        exit(1)

    for f in files:
        handle_file(f, options)

if __name__ == '__main__':
    main()

