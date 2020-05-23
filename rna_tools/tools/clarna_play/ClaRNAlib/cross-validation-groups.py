#!/usr/bin/env python
import numpy as np
import random
from utils import *

INP_FN = FileNamesObject.groups_fn(setname="training",reduced=True)
GROUPS_NUM = 5

random.seed(12345)

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
    print n_type, key, len(elems)
    random.shuffle(elems)
    
    s = np.array_split(np.array(elems), GROUPS_NUM)
    for i,sub_elems in enumerate(s):
        groups[i][key] += sub_elems.tolist()
    for i,sub_elems in enumerate(s):
        groups[i][rev_key] += [DoubletDescTool.reverse_d_id(x) for x in sub_elems]

# make all sets unique
for i in range(GROUPS_NUM):
    for key,elems in groups[i].items():
        old_count = len(elems)
        groups[i][key] = list(set(elems))
        new_count = len(groups[i][key])
        print i,key,new_count

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
        fn = FileNamesObject.groups_fn(setname=setname,reduced=True,cross_validation_num=i+1)
        save_json(fn, res_data)