#!/usr/bin/env python
import numpy as np
import random
from utils import *

DATA_FROM_GCH = [[9234, ["3u5h"]], [8443, ["3r8t"]], [8368, ["1s72"]], [8268, ["2j03"]], [7531, ["2zjr"]], [9271, ["3u5f", "1fjg"]], [8899, ["2nz4", "3g78", "2qbz", "1gid", "3kfu", "3rkf", "3akz", "3hhn", "3v7e", "3r1c"]], [9019, ["1m5o", "3p22", "1lng", "3dir", "3owz", "3rg5", "3add", "3d0u", "3adc", "1i9v", "1g1x", "1et4", "3cun", "2bte", "1u9s", "3q3z", "3d2v", "2hw8", "2gdi"]], [8952, ["2qwy", "2cv1", "2oiu", "2ho7", "3npn", "3ivn", "2xd0", "2zm5", "3fo4", "2csx", "3slq", "2quw", "1gax", "2azx", "2zzn", "3ol9", "3ol7", "3ol8", "1qa6", "2fmt", "3eph", "1il2", "3rw6", "1efw"]], [8781, ["1l9a", "3f4h", "2qux", "1zx7", "3sux", "1qcu", "3gx5", "3sd3", "3mut", "3fu2", "1s03", "2vpl", "1yyw", "2y9h", "1nuv", "1duq", "1xpe", "3bnq", "1q96", "3am1", "3a6p", "3bnp", "2nok", "1i6u", "2pn4", "3p59", "1h3e", "3la5", "2hoj", "2zuf", "1qrs", "1d4r", "1yls", "1feu", "1qu2", "1u0b", "1y26"]], [8895, ["1j1u", "1kxk", "3nvi", "3l0u", "3hax", "1qbp", "2zzm", "1vby", "3amt", "1mji", "3ova", "2ez6", "2tra", "3dh3", "1qf6", "1tfw", "2g3s", "3ouy", "3hjw", "1t0d", "2du3", "2i82", "3b5f", "1h4s", "3egz", "1jbr", "3moj", "1ser", "2p7d", "2d2l", "2gjw", "1b23", "1yfg", "2oeu", "1qc0", "1mzp", "3snp", "2ozb", "364d", "1f7v", "1f7y", "2pxv", "2dr2", "2rfk", "1sdr", "2dlc", "2ply", "3nj7", "1nlc", "3ftf", "2hvy", "3e5c", "3siu", "2nug", "2nue"]], [8812, ["1vfg", "1ooa", "2fk6", "157d", "1xjr", "2qek", "3ftm", "280d", "2o3v", "3bnl", "2fqn", "2zko", "1mwl", "2f8s", "2o3x", "3b31", "1z7f", "1r3e", "2bh2", "3mei", "1rpu", "299d", "1o9m", "3loa", "3ks8", "2q1o", "2zi0", "3r1d", "1r9f", "1yz9", "429d", "1urn", "1yzd", "3sj2", "1k8w", "2zy6", "1ntb", "3q50", "1kh6", "1dfu", "2bq5", "1di2", "1t0e", "2dr8", "1sa9", "3s49", "2oe5", "3bt7", "1a9n", "405d", "1rlg", "1zev", "420d", "3iab", "2b3j", "1rna", "2jlt", "1f1t", "406d", "1l2x", "2ao5", "1f27", "433d", "2bu1", "2ec0", "2awe", "3dvz", "1q2r", "2r22", "1ec6", "3og8", "1xok", "1e7k", "1duh", "361d", "1csl", "1zbh", "2az2", "1jid", "3lrn", "3r9x", "2i91", "5msf", "1saq", "2a43", "1u1y", "2w89", "466d", "3glp", "2e9t", "3cgp", "3tmi", "2val", "1i9x", "2xlj", "2bgg", "3lrr"]], [4409, ["353d", "1mhk", "435d", "2xsl", "3ncu", "2pjp", "3cgs", "1j9h", "397d", "1jzv", "3cgr", "2r20", "3ptx", "1lnt", "409d", "1kd5", "2iz9", "1fuf", "2c4z", "2izm", "3fs0", "2ann", "3ibk", "3cjz", "3oin", "422d", "1dqh", "7msf", "1a34", "2c50", "3r9w", "2gic", "3r2d", "1sds", "3bso", "377d", "3l25", "259d", "1gtf", "2izn", "1f8v", "1zse", "6msf", "472d", "354d", "3eqt", "3bsn", "3qrp", "402d", "1j8g", "438d", "2y8y", "2e9z", "2f8k", "439d", "3jxr", "2xli", "2v7r", "3bnt", "2y8w", "3gvn", "2vuq", "1ytu", "2v6w", "1cvj", "3mqk", "1j6s", "2g91", "3r1e", "3nd3", "1zh5", "2e9r", "2gxb", "3czw", "2db3", "1n35", "2a0p", "1yvp", "1kfo", "1zdj", "1m8v", "1utd", "3kmq", "1gtn", "1utf", "3koa", "413d", "3o7v", "3nl0", "1l3z", "2q1r", "3kms", "255d", "1wne", "283d", "3mij", "3m85", "2jlx", "2b2d", "1uvm", "2jlw", "3m7n", "3nj6", "2bx2", "2xgj", "2jlz", "2jlu", "1bmv", "2grb", "2j0s", "3mj0", "1si3", "2vnu", "3fht", "1m8y", "2atw", "3l26", "2ix1", "2r7w", "3p6y", "2von", "2vod", "1b7f", "3boy", "1pgl", "2bbv", "1n38", "3iev", "3rc8", "1osu", "1rxa", "1knz", "1g2j", "1wmq", "3sqw", "1p79", "2r7t", "3i5x", "3pey", "3mdi", "2xs5", "2xs2", "1uvj", "1wpu", "3hga", "3aev", "1av6", "2po1", "3qgc", "1fxl", "3nma", "2asb", "1g2e", "1mdg", "3mdg", "2r7v", "1n1h", "3k5y", "3k5z", "3k5q", "3k64", "3k61", "3k62", "333d", "2jea", "2xzo", "3nnh", "3ice", "3pf5", "3pf4", "3i5y", "3g9y", "3bx2", "1ddl", "3r2c", "2xs7", "1uvi", "1uvn", "3q0s", "3rer", "3ie1", "3o8c", "3o8r", "3d2s", "2vop", "3q0q", "3hsb", "3bsb", "3bsx", "3ahu", "1nb7", "2g4b", "3nna", "3nnc", "3u2e", "3er9", "1i5l", "2x1f", "1h2c", "1rmv", "2a8v", "3nmr", "3rtj", "2q66", "1m8w", "2tmv", "3k49", "3gib", "2xnr", "1uvl", "3q0o", "2a1r", "1pvo", "3q0p", "3q0r", "3gpq", "2voo", "2fz2", "3bx3", "1kq2"]]]

# self-check
pdbids_in_training = set(open("data-training.txt").read().split("\n"))
p2g = {}
for group_num,(count,elems) in enumerate(DATA_FROM_GCH):
    for pdbid in elems:
        if not pdbid in pdbids_in_training:
            print "file %s not in training" % pdbid
            continue
        assert pdbid in pdbids_in_training, "file %s not in training set" % pdbid
        if p2g.has_key(pdbid):
            print "file %s already in group %s" % (pdbid,p2g[pdbid])
            continue
        assert not p2g.has_key(pdbid)
        p2g[pdbid]=group_num

INP_FN = FileNamesObject.groups_fn(setname="training",reduced=True)
GROUPS_NUM = len(DATA_FROM_GCH)

def filter_doublets(all_doublets, f_doublets):
    return [x for x in all_doublets if x in f_doublets]

data = load_json(INP_FN)

groups = [{} for i in range(GROUPS_NUM)]
for n_type in N_TYPES:
    key = "all/all/%s"%n_type
    rev_key = "all/all/%s"%DoubletDescTool.reverse_n_type(n_type)

    for i in xrange(GROUPS_NUM):
        for k in (key,rev_key):
            if not groups[i].has_key(k):
                groups[i][k] = []

    elems = data.get(key,[])
    s = [[] for i in xrange(GROUPS_NUM)]
    for d_id in elems:
        pdbid = d_id.split(":")[0].lower()
        if p2g.has_key(pdbid):
            s[p2g[pdbid]].append(d_id)

    for i,sub_elems in enumerate(s):
        groups[i][key] += sub_elems
    for i,sub_elems in enumerate(s):
        groups[i][rev_key] += [DoubletDescTool.reverse_d_id(x) for x in sub_elems]


for i in range(GROUPS_NUM):
    training = {}
    validate = {}
    for n_type in N_TYPES:
        key = "all/all/%s"%n_type
        training[n_type]=[]
        validate[n_type]=[]
        for j in range(GROUPS_NUM):
            if i!=j:
                training[n_type] += groups[j][key]
            else:
                validate[n_type] += groups[j][key]
        training[n_type]=set(training[n_type])
        validate[n_type]=set(validate[n_type])
        # exclude ZZ files from validate set
        zz_elems = set([x for x in validate[n_type] if re.match('^ZZ',x)])
        validate[n_type] = validate[n_type].difference(zz_elems)
        print "i=%d n_type=%s training=%d validate=%d removed zz_elems=%d" % (i,n_type,len(training[n_type]),len(validate[n_type]),len(zz_elems))
    t_data = {}
    v_data = {}
    
    for key,elems in data.items():
        n_type = key.split("/")[-1]
        t_data[key] = filter_doublets(elems, training[n_type])
        v_data[key] = filter_doublets(elems, validate[n_type])
        print "i=%d key=%s elems=%d t_data=%d v_data=%d" % (i,key,len(set(elems)),len(t_data[key]),len(v_data[key]))
        # assert "/other/" in key or len(elems)==(len(t_data[key])+len(v_data[key]))

    for setname,res_data in (("training",t_data),("bench",v_data)):
        fn = FileNamesObject.groups_fn(setname=setname,reduced=True,cross_validation_num=chr(ord('A')+i))
        save_json(fn, res_data)