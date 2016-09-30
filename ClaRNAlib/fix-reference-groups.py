#!/usr/bin/env python
#
# fixes to reference groups
#
import os
from utils import *
from distances import doublet_params_dict, expected_strand_orient

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""add new files to expert directory""")

    (options, args)  = parser.parse_args()
    return (parser, options, args)

class Fixer(object):

    def __init__(self):
        self.expert = {'not-ref':{},'fuzzy':{},'ref':{}}

    def load(self,comment):
        for d in 'not-ref','fuzzy','ref':
            for fn in os.listdir('expert/%s'%d):
                if re.match('^.*%s.json'%comment,fn):
                    data = load_json("expert/%s/%s"%(d,fn))
                    print fn, len(data), data[0:10]
                    self._add(d,fn.replace(".json",""),data)

    def _add(self,key1,key2,doublets):
        if not self.expert[key1].has_key(key2):
            self.expert[key1][key2] = set()
        for d_id in doublets:
            self.expert[key1][key2].add(d_id)

    def add_to_ref(self,sc,desc,n_type,doublets,comment="",add_symmetric=True):
        self._add('ref', "%s_%s_%s_%s" % (sc,n_type,desc,comment), doublets)
        if add_symmetric:
            rev_n_type = DoubletDescTool.reverse_n_type(n_type)
            rev_desc = DoubletDescTool.reverse_desc(desc)
            rev_doublets = [DoubletDescTool.reverse_d_id(x) for x in doublets]
            self._add('ref', "%s_%s_%s_%s" % (sc,rev_n_type,rev_desc,comment), rev_doublets)


    def from_ref_to_fuzzy(self,sc,desc,n_type,doublets,comment="",add_symmetric=True):
        self._add('not-ref', "%s_%s_%s_%s" % (sc,n_type,desc,comment), doublets)
        self._add('fuzzy', "%s_%s_%s_%s" % (sc,n_type,desc,comment), doublets)
        if add_symmetric:
            rev_n_type = DoubletDescTool.reverse_n_type(n_type)
            rev_desc = DoubletDescTool.reverse_desc(desc)
            rev_doublets = [DoubletDescTool.reverse_d_id(x) for x in doublets]
            self._add('not-ref', "%s_%s_%s_%s" % (sc,rev_n_type,rev_desc,comment), rev_doublets)
            self._add('fuzzy', "%s_%s_%s_%s" % (sc,rev_n_type,rev_desc,comment), rev_doublets)

    def save(self):
        for key1 in ('ref','not-ref','fuzzy'):
            for key2 in sorted(self.expert[key1].keys()):
                doublets = sorted(list(self.expert[key1][key2]))
                print key1,key2,doublets
                save_json("expert/%s/%s.json"%(key1,key2), doublets, indent=2)

def tmp_gen_script():
    eval_fn = FileNamesObject.eval_fn(setname="bench",t="eval1")
    json = load_json(eval_fn)
    # move from ref to fuzzy
    for group,desc,n_type in (
        ('ref-ok','HH_cis','CG'),
        ('ref-ok-fuzzy','HH_cis','CG'),
        #
        ('ref-undetected','HW_cis','UC'),
        ('ref-undetected','HW_tran','AG'),
        #
        ('ref-undetected','SW_cis','GC'),
        ('ref-ok-fuzzy','SW_cis','GC'),
        #
        ('ref-undetected','HW_tran','UA'),
        ('ref-ok-fuzzy','HW_tran','UA'),
        #
        ('ref-undetected','SH_tran','GC'),
        ('ref-ok-fuzzy','SH_tran','GC'),
        #
        ('ref-ok-fuzzy','SW_cis','UU'),
        #
        ('ref-ok','WH_cis','AC'),
        ('ref-undetected','WH_cis','AC'),
        ('ref-ok-fuzzy','WH_cis','AC'),
        #
        ('ref-undetected','WS_tran','AC'),
        ('ref-ok','WS_tran','AC'),
        #
        ('ref-undetected','WW_tran','UG'),
        ('ref-ok-fuzzy','WW_tran','UG'),
    ):
        key = "evaluation/%(group)s/%(desc)s/%(n_type)s" % locals()
        doublets = json[key]
        print '    f.from_ref_to_fuzzy("bp","%(desc)s","%(n_type)s",%(doublets)s,"bad_match")' % locals()

def tmp_gen_from_eval():

    def compute_doublet_params(ids, params_type):
        dp = __import__("doublet-params")
        
        R_ATOMS = ['N1','N2','N3','N4','N6','N7','P']
        R_ATOMS += ['C2','C4','C5','C6','C8',"C1'"]
        R_ATOMS += ['O2']
        R_ATOMS += ["O2'","O3'","O4'"]
        R_ATOMS += ["OP1","OP2","O5'","NEXT:O3'"]
        R_ATOMS += ['N1','C6','O6','C5','C4','N3','C2','N2','N7','C8','N9']

    
        dd = DoubletsDict("gc-data",reduced_atoms=R_ATOMS)
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

    f = Fixer()
    f.load("bad_params")
    f.load("wrong_orient")
    
    for fn in [FileNamesObject.eval_fn(setname="training",t="eval1"),FileNamesObject.eval_fn(setname="bench",t="eval1")]:
        json = load_json(fn)
        json = dict([(k,v) for k,v in json.items() if re.match('^evaluation/(ref-undetected|ref-diff)/.._(cis|tran)/.*',k)])
        ids = sum(json.values(),[])
        print len(ids), ids[0:100]
        
        all_params = compute_doublet_params(ids,'bp')
        
        limits = {'nn_ang_norm': {'min':0.0,'max':65.0}, 
                  'dist': {'min':4.0,'max':8.5}, 
                  'n1cc_ang': {'min': 50.0,'max':140.0}, 
                  'n2cc_ang': {'min': 50.0,'max':140.0}
                  }
        
        ids2 = filter_by_params(ids, limits, all_params) 
        bad_doublets = set(ids).difference(set(ids2))
        
        print "bad doublets: %d" % len(bad_doublets)
        for k,v in json.items():
            n_type = k.split("/")[-1]
            desc = k.split("/")[-2]
            b = bad_doublets.intersection(set(v))
            if len(b)>0:
                print desc, n_type, len(b), "bad_params"
                f.from_ref_to_fuzzy("bp",desc,n_type,b,"bad_params")
            good_strand_orient = expected_strand_orient(desc)
            bad_orient = []
            for d_id in v:
                if all_params[d_id]['strand_orient'] != good_strand_orient:
                    bad_orient.append(d_id)
            if len(bad_orient)>0:
                print desc, n_type, len(bad_orient), "wrong_orient"
                f.from_ref_to_fuzzy("bp",desc,n_type,bad_orient,"wrong_orient")
    f.save()

def main():
    (parser,options,_args) = parse_args()

    f = Fixer()
    f.from_ref_to_fuzzy("bp","HH_cis","CG",['2KRL:A25:A28'],"bad_match")
    f.from_ref_to_fuzzy("bp","HH_cis","CG",['2HGP:D16:D18','3E1B:B1134:B1137'],"bad_match")
    f.from_ref_to_fuzzy("bp","HH_cis","CG",['3J14:A87:B951'],"bad_match")
    f.from_ref_to_fuzzy("bp","HW_cis","UC",['3HUZ:A2027:A2080'],"bad_match")
    f.from_ref_to_fuzzy("bp","HW_cis","UG",['2HGQ:A1113:A1094','2WH1:W58:W56'],"bad_match")
    f.from_ref_to_fuzzy("bp","HW_tran","AG",['2B66:A305:A383', '2B9P:A305:A383', '2QBE:B1028:B1090', '2Y0V:A1569:A1582', '2Y0X:A1569:A1582', '2YKR:A1439:A1459'],"bad_match")
    f.from_ref_to_fuzzy("bp","HW_tran","AG",['1VSA:w1747:w1743','2HGJ:A1744:A1746','2HGU:A1810:A1825',
          '2OW8:y1057:y1054','3CXC:0293:0333','3E1B:B2452:B2454','3E1D:B250:B246'],"bad_match")
    f.from_ref_to_fuzzy("bp","SW_cis","GC",['3I8H:A1104:A1124', '3V26:A1102:A1122', '3V28:A1098:A1118'],"bad_match")
    f.from_ref_to_fuzzy("bp","SW_cis","GC",['1Z58:21257:21270', '2F4V:A975:A1016', '2UXB:A1104:A1124', '2WRN:A1025:A1190', '2WRQ:A1025:A1190', '2X9R:A1104:A1124', '2X9T:A1104:A1124', '2XG1:W5:W68', '2Y0U:A976:A1016', '2Y0W:A976:A1016', '2Y0Y:A976:A1015', '2Y12:A976:A1015', '3HUW:A968:A1008', '3HUY:A1104:A1124', '3HUY:A1025:A1190', '3KCR:81055:81103', '3KIU:a976:a1015', '3KIX:a976:a1015', '3ZVP:A1133:A1140', '4DR3:A1104:A1124', '4GAQ:A1046:A1212'],"bad_match")
    f.from_ref_to_fuzzy("bp","HW_tran","UA",['3UZ6:A805:A844', '3UZL:A805:A844'],"bad_match")
    f.from_ref_to_fuzzy("bp","HW_tran","UA",['3OFQ:A1418:A1532', '3TVF:A805:A844'],"bad_match")
    f.from_ref_to_fuzzy("bp","SH_tran","GC",['2WDM:Y26:Y42'],"bad_match")
    f.from_ref_to_fuzzy("bp","SH_tran","GC",['2WRO:A1574:A1577', '2WRR:A1574:A1577', '2XFZ:A1106:A1126', '3D5B:A1552:A1555', '3D5D:A1552:A1555'],"bad_match")
    f.from_ref_to_fuzzy("bp","SH_tran","GC",['1VQ4:02365:02367','1VQ6:02365:02367','1YIT:02365:02367',
          '2OTL:02365:02367','3CC7:02365:02367','3G71:02365:02367'],"bad_match")
    f.from_ref_to_fuzzy("bp","SW_cis","UU",['1J5A:A2183:A2190', '1JZX:A2183:A2190', '1JZY:A2183:A2190', '1JZZ:A2183:A2190', '1K01:A2183:A2190', '2I2T:B1153:B1156', '2KUV:A20:A25', '2KUW:A20:A25', '2QAO:B1153:B1156', '2QBC:B1153:B1156', '3J0T:B2806:B2890', '3J14:B2806:B2890', '4A18:1451:1525', '4A19:1451:1525', '4A1B:1451:1525', '4A1D:1451:1525'],"bad_match")
    f.from_ref_to_fuzzy("bp","WH_cis","AC",['2WW9:D13:D21'],"bad_match")
    f.from_ref_to_fuzzy("bp","WH_cis","AC",['2WWB:D13:D21'],"bad_match")
    f.from_ref_to_fuzzy("bp","WH_cis","AC",['3IZT:B2198:B2225'],"bad_match")
    f.from_ref_to_fuzzy("bp","WH_cis","UC",['3I1P:A2456:A2458'],"bad_match")
    f.from_ref_to_fuzzy("bp","WS_tran","AC",['3UYE:A2148:A2185', '3UZ1:A2148:A2185','1S1I:31116:32451'],"bad_match")
    f.from_ref_to_fuzzy("bp","WS_tran","AC",['2QOV:B1113:B2005', '2QOX:B1113:B2005'],"bad_match")
    f.from_ref_to_fuzzy("bp","WS_tran","AC",['1VOR:B2217:B2213','1VOU:B2217:B2213','2HOM:A49:A65',
          '3IZV:A693:A795','3IZV:A899:A809','3J0V:A899:A809','3J0X:A1200:A961','3J0X:A899:A809',
          '3UXR:B52:B30','3V27:B52:B30'],"bad_match")
    f.from_ref_to_fuzzy("bp","WW_tran","UG",['1Z31:A14:A17', '2AVY:A1130:A1133', '3FIK:B199:B249', '4A2I:A1130:A1133', '4ADV:A1130:A1133'],"bad_match")
    f.from_ref_to_fuzzy("bp","WW_tran","UG",['2XZM:A1149:A1163', '2XZN:A1149:A1163', '2ZJP:X1305:X1309', '3HUW:A1416:A1419', '3IZT:B199:B249', '3IZU:B199:B249', '3J0L:g16:g30', '3J0O:g16:g30', '3J0P:g16:g30', '3J0Q:g16:g30', '3JYX:5695:5700', '3ZVP:A1942:A1945'],"bad_match")
    f.from_ref_to_fuzzy("bp","WW_tran","UG",['1I97:A1115:A1118','1VQ8:02486:02381','1ZJW:B52:B15',
          '3E1A:11130:11133','3FIH:V55:V18','3J0U:A1449:A1452','3J0V:A1449:A1452',
          '3J0X:A1449:A1452'],"bad_match")
    f.from_ref_to_fuzzy("bp","WW_tran","GU",['3FIH:V18:V55'],"bad_match")
    
    f.from_ref_to_fuzzy("bp","WW_tran","UG",["1M82:A10:A13", "2OGM:01633:01636", "2X9R:A1425:A1428", "3I1S:A338:A341", "3J10:A1449:A1452", 
    "3J11:B1691:B1694", "3KIQ:a1425:a1428", "3KIX:a1425:a1428", "4E8R:A341:A344"],"bad_match")
    
    f.from_ref_to_fuzzy("bp","WH_cis","CA",["1VOW:A41:A43", "1VOY:A41:A43", "1VP0:A41:A43", "2O43:A57:A72","3IZT:B1294:B1274"],"bad_match")
    f.from_ref_to_fuzzy("bp","WH_tran","AG",["3IZT:B1663:B1991", "3JQ4:A466:A25", "3UYE:A1149:A1102", "3UYF:A260:A256", "3UZ1:A1149:A1102"])
    
    f.from_ref_to_fuzzy("bp","HW_cis","AC",["1Z58:272:257", "2GYC:0434:0432", "1VOR:A43:A41", "1VOU:A43:A41"],"bad_match")
    f.from_ref_to_fuzzy("bp","SW_tran","CA",["1VOW:B2213:B2217", "1VOY:B2213:B2217", "1VP0:B2213:B2217", "3I9E:A2056:A2563", 
          "3J10:A809:A899", "4A1B:1901:1906", "3V29:B30:B52"],"bad_match")
    f.from_ref_to_fuzzy("bp","HW_tran","CU",["3FWO:A2096:A1997"],"bad_match")
    f.from_ref_to_fuzzy("bp","WH_tran","CU",["2UUA:A342:A51","3U5B:2423:251"],"bad_match")
    f.from_ref_to_fuzzy("bp","WH_tran","CU",["2GYC:0300:0239", "2QBE:B2359:B2330", "2QBG:B2359:B2330", "2QBH:A347:A51", 
          "2QBJ:A347:A51", "2WWQ:B2555:B1954", "3I1S:A347:A51", "3I1Z:A350:A54", "3I21:A347:A51", "3JYX:52013:52003",
          "1I95:A345:A54", "3UZ1:A2889:A2846", "3UZ4:A342:A51"],"bad_match")
    f.from_ref_to_fuzzy("bp","WS_cis","UU",["2GYA:01375:01381", "2GYC:01375:01381", "3JYX:5425:5621"],"bad_match")
    f.from_ref_to_fuzzy("bp","HS_cis","UU",["3IZE:A884:A883", "4A1B:12627:12623"],"bad_match")

    f.from_ref_to_fuzzy("bp","HS_cis","CU",['2ZJR:X2308:X2306', '3E1D:B435:B433', '3KNI:A2609:A2607'],"border_match") 
    f.from_ref_to_fuzzy("bp","SS_cis","AG",['3OHZ:A358:A72', '3OI1:A358:A72'],"border_match") 
    f.from_ref_to_fuzzy("bp","SS_cis","CG",['3J0Y:A25:A115'],"border_match") 
    f.from_ref_to_fuzzy("bp","SS_cis","GA",['3OHZ:A72:A358', '3OI1:A72:A358'],"border_match") 
    f.from_ref_to_fuzzy("bp","SW_cis","GA",['2WH1:A1217:A1214', '3BBO:A171:A1330', '3J0Y:B951:B2266', '3J0Z:A1240:A1237', '3KNI:A72:A358', '3ZVO:A1217:A1214'],"border_match") 
    f.from_ref_to_fuzzy("bp","WS_cis","AG",['2WH1:A1214:A1217', '3IZW:A900:A768', '3J0Y:B2266:B951', '3J0Z:A1237:A1240', '3ZVO:A1214:A1217'],"border_match") 
    f.from_ref_to_fuzzy("bp","WS_cis","AU",['2VHN:B1053:B1045'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","AC",['2E5L:A1472:A1381', '2J01:A1705:A695', '2J03:A1705:A695', '3I8I:A1811:A734', '3R8S:B108:B10'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","AU",['1VOR:B1135:B1125', '1VOU:B1135:B1125'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","CC",['1VOR:A17:A63', '1VOU:A17:A63'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","CG",['2J00:A1002:A1009', '2J02:A1002:A1009', '2WH2:A2860:A2834'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","CU",['2WH1:A1117:A1112', '3O58:318:1402'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","GC",['2J01:A2720:A2746', '2J03:A2720:A2746', '2WH2:A2834:A2860'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","GU",['1I97:A869:A880', '3IZV:A1418:A1480'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","UA",['3E1B:B1647:B2008'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","UG",['1I97:A880:A869', '2OGN:02068:02002', '3IZV:A1480:A1418'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_tran","CG",['3E1A:E47:E14'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_tran","GC",['3E1A:E14:E47'],"border_match") 
    f.from_ref_to_fuzzy("bp","HS_cis","AC",['3UYD:A364:A362'],"border_match") 
    f.from_ref_to_fuzzy("bp","HS_cis","CU",['3I8F:A2530:A2528'],"border_match") 
    f.from_ref_to_fuzzy("bp","HS_cis","GG",['2QB9:A249:A246', '2QBB:A249:A246'],"border_match") 
    f.from_ref_to_fuzzy("bp","HW_cis","CC",['3EOG:A357:A374', '3J14:B866:B864'],"border_match") 
    f.from_ref_to_fuzzy("bp","HW_cis","CU",['3DF2:B163:B161', '3DF4:B163:B161'],"border_match") 
    f.from_ref_to_fuzzy("bp","HW_cis","GG",['2D17:A11:A9'],"border_match") 
    f.from_ref_to_fuzzy("bp","SH_cis","CA",['3UYD:A362:A364'],"border_match") 
    f.from_ref_to_fuzzy("bp","SH_cis","GG",['2QB9:A246:A249', '2QBB:A246:A249'],"border_match") 
    f.from_ref_to_fuzzy("bp","SH_cis","UC",['3I8F:A2528:A2530'],"border_match") 
    f.from_ref_to_fuzzy("bp","SH_tran","UC",['1YKV:B30:B13'],"border_match") 
    f.from_ref_to_fuzzy("bp","SS_cis","AG",['1OND:0171:0226', '2X9S:A1405:A170', '2X9U:A1405:A170', '3OI3:A1314:A168', '3UZ8:A1406:A171', '3UZH:A365:A77'],"border_match") 
    f.from_ref_to_fuzzy("bp","SS_cis","GA",['1OND:0226:0171', '3OI5:A168:A1314', '4ABS:A72:A364'],"border_match") 
    f.from_ref_to_fuzzy("bp","SS_cis","UG",['3MR8:A1341:A1209', '3MS0:A1341:A1209'],"border_match") 
    f.from_ref_to_fuzzy("bp","SW_cis","CA",['1HNZ:A1380:A1471'],"border_match") 
    f.from_ref_to_fuzzy("bp","SW_cis","GA",['1VS8:B186:B1345', '2AWB:B186:B1345', '2B9P:A77:A370', '2I2T:B186:B1345', '2I2V:B186:B1345', '2J28:B186:B1345', '2QP0:A1236:A1233', '2WH3:A1217:A1214', '2WWQ:B186:B1364', '2Y15:A72:A364', '3I1Q:A1239:A1236', '3I1Q:A1258:A1259', '3IYQ:A211:A216', '3IZT:B2591:B1965', '3J01:8186:81364', '3J0W:B69:B70', '3J11:B951:B2266', '3J12:B951:B2266', '3KNJ:A1208:A1205', '3KNM:A558:A1167', '3OFR:A186:A1358', '3SGF:A186:A1357', '4ABR:A1217:A1214'],"border_match") 
    f.from_ref_to_fuzzy("bp","SW_cis","GC",['2X9R:A1011:A1000', '2X9T:A1011:A1000', '4DHA:A279:A276'],"border_match") 
    f.from_ref_to_fuzzy("bp","SW_tran","UA",['1YL4:A933:A1177'],"border_match") 
    f.from_ref_to_fuzzy("bp","WH_cis","CC",['3EOG:A374:A357'],"border_match") 
    f.from_ref_to_fuzzy("bp","WH_cis","UC",['2X9S:A2525:A2527'],"border_match") 
    f.from_ref_to_fuzzy("bp","WS_cis","AG",['1VS8:B1345:B186', '2AWB:B1345:B186', '2I2T:B1345:B186', '2I2V:B1345:B186', '2J28:B1345:B186', '2QP0:A1233:A1236', '2WH3:A1214:A1217', '2Y15:A364:A72', '3I1Q:A1236:A1239', '3IZT:B1965:B2591', '3J01:81364:8186', '3KNJ:A1205:A1208', '3KNM:A1167:A558', '3OFR:A1358:A186', '3SGF:A1357:A186', '4ABR:A1214:A1217', '4FAX:A261:A7'],"border_match") 
    f.from_ref_to_fuzzy("bp","WS_tran","AU",['1YL4:A1177:A933', '2WDG:Y20:Y7'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","AC",['1PNU:9108:99', '1PNY:9108:99', '2XKV:B68:B27', '3FIN:A1757:A695', '3I1O:A200:A209', '3J19:B108:B10', '3OAR:A200:A209', '3ORA:A1494:A1398', '4DR5:A1471:A1380'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","AU",['1VOW:B1135:B1125', '1VOY:B1135:B1125', '1VP0:B1135:B1125'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","CC",['1VOW:A17:A63', '1VOY:A17:A63', '1VP0:A17:A63', '2HGI:A574:A623'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","CG",['2WH4:A2860:A2834', '2WWQ:B2176:B2120', '3FIC:A1002:A1009', '3J0T:B2829:B2815', '3KNM:A2764:A2738'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","CU",['2WH3:A1117:A1112', '3O5H:318:1402', '3OI2:A1117:A1112', '3OI4:A1117:A1112'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","GC",['2WH4:A2834:A2860', '2WWQ:B2120:B2176', '3FIN:A2800:A2826', '3KIR:A1099:A1082', '3KIT:A1099:A1082', '3KNM:A2738:A2764', '3MRZ:A2827:A2853', '3MS1:A2827:A2853', '4ABR:W42:W28'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","GG",['3HUW:A1100:A1112'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","GU",['1Q7Y:A61:A97', '2VHP:A645:A584', '3I55:0632:0890', '3J0W:B2233:B2084', '3V28:A1070:A1057'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","UA",['4AQY:A1181:A1032', '4DR6:A1181:A1032'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","UC",['3OI2:A1112:A1117', '3OI4:A1112:A1117'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","UG",['3I55:0890:0632'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_cis","UU",['3KCR:82488:82490'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_tran","CG",['4DR5:A182:A118'],"border_match") 
    f.from_ref_to_fuzzy("bp","WW_tran","GC",['2NQP:F6:F21', '3UYD:A256:A261', '3UZG:A256:A261', '4DR5:A118:A182'],"border_match") 
    f.from_ref_to_fuzzy("bp","SH_cis","GG",["2QBH:A246:A249", "2QBJ:A246:A249", "3DF1:A246:A249", "3DF3:A246:A249"],"border_match")


    f.add_to_ref("bp","HH_tran","CC",['ZZ01:A2:B2'],"from_JS") 
    f.add_to_ref("bp","SH_cis","GC",['ZZ02:A1:B1'],"from_JS") 
    f.add_to_ref("bp","SH_cis","AG",['ZZ03:A1:B1'],"from_JS") 
    f.add_to_ref("bp","HS_tran","GA",['ZZ04:A2:B2'],"from_JS") 
    f.add_to_ref("bp","HS_tran","GC",['ZZ05:A2:B2'],"from_JS") 
    f.add_to_ref("bp","HH_tran","GG",['ZZ06:A2:B2'],"from_JS") 
    f.add_to_ref("bp","HS_tran","GU",['ZZ07:A2:B2'],"from_JS") 
    f.add_to_ref("bp","SH_cis","AU",['ZZ08:A1:B1'],"from_JS") 
    f.add_to_ref("bp","HS_tran","UA",['ZZ09:A2:B2'],"from_JS") 
    f.add_to_ref("bp","HS_tran","UC",['ZZ10:A2:B2'],"from_JS") 
    f.add_to_ref("bp","HS_tran","UU",['ZZ11:A2:B2'],"from_JS")
    f.add_to_ref("bp","HW_cis","UG",['ZZ12:A0:B1'],"morphed_from_UA")
    f.add_to_ref("bp","WH_cis","CA",['ZZ13:A0:B1'],"morphed_from_CG")

    f.add_to_ref("bp","HH_cis","AA",['ZZ20:A0:B1', 'ZZ21:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HH_cis","AC",['ZZ22:A0:B1', 'ZZ23:A0:B1', 'ZZ24:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HH_cis","AU",['ZZ25:A0:B1', 'ZZ26:A0:B1', 'ZZ27:A0:B1', 'ZZ28:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HH_cis","AG",['ZZ29:A0:B1', 'ZZ2A:A0:B1', 'ZZ2B:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HH_cis","CG",['ZZ2C:A0:B1', 'ZZ2D:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HS_tran","CG",['ZZ2E:A0:B1', 'ZZ2F:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HW_cis","AA",['ZZ30:A0:B1', 'ZZ31:A0:B1', 'ZZ32:A0:B1', 'ZZ33:A0:B1', 'ZZ34:A0:B1', 'ZZ35:A0:B1', 'ZZ36:A0:B1', 'ZZ37:A0:B1', 'ZZ38:A0:B1', 'ZZ39:A0:B1', 'ZZ3A:A0:B1', 'ZZ3B:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HW_tran","CU",['ZZ3C:A0:B1', 'ZZ3D:A0:B1', 'ZZ3E:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","HW_tran","CG",['ZZ3F:A0:B1', 'ZZ40:A0:B1', 'ZZ41:A0:B1', 'ZZ42:A0:B1', 'ZZ43:A0:B1', 'ZZ44:A0:B1', 'ZZ45:A0:B1', 'ZZ46:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","SH_cis","CU",['ZZ47:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","SS_cis","AA",['ZZ48:A0:B1', 'ZZ49:A0:B1', 'ZZ4A:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","SS_cis","AC",['ZZ4B:A0:B1', 'ZZ4C:A0:B1', 'ZZ4D:A0:B1', 'ZZ4E:A0:B1', 'ZZ4F:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","SS_cis","AU",['ZZ50:A0:B1', 'ZZ51:A0:B1', 'ZZ52:A0:B1', 'ZZ53:A0:B1', 'ZZ54:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","SS_tran","AA",['ZZ55:A0:B1', 'ZZ56:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","SW_tran","CU",['ZZ57:A0:B1', 'ZZ58:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","WH_cis","CU",['ZZ59:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","WH_tran","AU",['ZZ5A:A0:B1', 'ZZ5B:A0:B1', 'ZZ5C:A0:B1', 'ZZ5D:A0:B1', 'ZZ5E:A0:B1', 'ZZ5F:A0:B1', 'ZZ60:A0:B1', 'ZZ61:A0:B1'],"morphed_from_normal")
    f.add_to_ref("bp","WS_tran","CU",['ZZ62:A0:B1'],"morphed_from_normal")

    # 2013-06-03
    f.add_to_ref("bp","HH_cis","CC",['ZZ70:A0:B1'],"morphed")
    f.add_to_ref("bp","HH_cis","CU",['ZZ71:A0:B1'],"morphed")
    f.add_to_ref("bp","HH_cis","UU",['ZZ72:A0:B1'],"morphed")
    f.add_to_ref("bp","HH_cis","GU",['ZZ73:A0:B1'],"morphed")
    f.add_to_ref("bp","HH_cis","GG",['ZZ74:A0:B1'],"morphed")
    f.add_to_ref("bp","HH_tran","UU",['ZZ75:A0:B1'],"morphed")
    f.add_to_ref("bp","HH_tran","GU",['ZZ76:A0:B1'],"morphed")
    f.add_to_ref("bp","HS_cis","GC",['ZZ77:A0:B1'],"morphed")
    f.add_to_ref("bp","HS_cis","GU",['ZZ78:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_cis","CC",['ZZ79:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_cis","CU",['ZZ7A:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_cis","UU",['ZZ7B:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_tran","AC",['ZZ7C:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_tran","AU",['ZZ7D:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_tran","CC",['ZZ7E:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_tran","CU",['ZZ7F:A0:B1'],"morphed")
    f.add_to_ref("bp","SS_tran","UU",['ZZ80:A0:B1'],"morphed")

    f.save()

if __name__=="__main__":
    main()
    # tmp_gen_script()
    # tmp_gen_from_eval()
