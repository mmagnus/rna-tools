#!/usr/bin/env python
import re
from utils import *
all = load_json("descriptions-dict.json")
mc = all['MC']
for k,v in sorted(mc.items()):
    if "stacking" in k:
        continue
    kk = k
    cl_code = v
    if "REV:" in kk:
        kk = k.replace("REV:","")
        cl_code = v[1]+v[0]+v[2:]
    orient = None
    if "cis" in kk:
        orient = "cis"
    elif "trans" in kk:
        orient = "tran"
    kk2 = kk[0]+kk[2]
    if orient is not None and re.match("^[HWS]{2}$",kk2):
        mc_code = kk2+"_"+orient
    else:
        mc_code = "UNK"
    if mc_code != cl_code:
       print "INVALID description! %s=%s should be %s" % (k,v,mc_code)
       if mc_code=='UNK':
          del all['MC'][k]
       else:
          all['MC'][k] = mc_code
save_json("descriptions-dict-new.json",all,indent=2)
