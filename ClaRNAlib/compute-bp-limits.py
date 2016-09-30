#!/usr/bin/env python
from utils import *
from distances import doublet_params_dict


R_ATOMS = ['P',"C1'","O2'","O3'","O4'"]
R_ATOMS += ["OP1","OP2","O5'","NEXT:O3'"]
# base atoms
R_ATOMS += ['N9', 'C4', 'N3', 'N1', 'C6', 'N6', 'C8', 'C5', 'C2', 'N7']
R_ATOMS += ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C6', 'C5']
R_ATOMS += ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C6', 'C5']
R_ATOMS += ['N9', 'C4', 'N3', 'N1', 'C6', 'O6', 'C8', 'C5', 'C2', 'N7', 'N2']


def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""some experiments with classifier generation""")
    parser.add_option("--groups", dest="groups",
                  help="groups filename", metavar="FILE", 
                  default=FileNamesObject.groups_fn(setname="training",reduced=True)+","+FileNamesObject.groups_fn(setname="bench",reduced=True))
    parser.add_option("--n-type", dest="n_type",
                  help="select n-type", metavar="NTYPE", default="AC")
    parser.add_option("--desc", dest="desc",
                  help="select classifier description", metavar="D", default="HW_tran")
    parser.add_option("--output-dir", dest="output_dir",
                  help="set output dir", metavar="DIR", default="/tmp/bp_limits")
    parser.add_option("--data-dir", dest="data_dir",
                  help="set data dir", metavar="DIR", default="gc-data")

    (options, args)  = parser.parse_args()
    return (parser, options, args)

def gen_limits(options):
    cc = __import__("compute-classifier")

    ref = set()
    for fn in options.groups.split(","):
        print >> sys.stderr, "processing %s" % fn
        json = load_json(fn)
        key1 = "classifier/bp/%s/%s" % (options.desc, options.n_type)
        key2 = "classifier/bp/%s/%s" % (DoubletDescTool.reverse_desc(options.desc), DoubletDescTool.reverse_n_type(options.n_type))
        ids1 = json.get(key1,[])
        ids2 = json.get(key2,[])
        ids2 = [DoubletDescTool.reverse_d_id(x) for x in ids2]
        ref.update(ids1)
        ref.update(ids2)

    if len(ref)==0:
        print "empty reference set!"
        return
    
    dd = DoubletsDict("gc-data",reduced_atoms=R_ATOMS)
    all_params = pickle_load(cc.fn("doublets-params",sc="bp",n_type=options.n_type))
    
    missing_ids = set(ref).difference(set(all_params.keys()))
    other_params = cc.compute_doublet_params(missing_ids, 'bp', options)
    for k,v in other_params.items():
        all_params[k] = v
    del other_params

    n = len(ref)
    limits = {}
    all_limits_str = ""
    for key in ['min_dist','dist','rot_ang','n2_z','nn_ang_norm']:
        values = [all_params[d_id][key] for d_id in ref]
        values.sort()
        shift_v = 0
        shift_start = None
        if key in ['rot_ang']:
            if min(values)<-20 and max(values)>200:
                v1 = [x for x in values if abs(-90-x)<abs(270-x)]
                v2 = [x for x in values if abs(-90-x)>abs(270-x)]
                assert len(v1)>0
                assert len(v2)>0
                shift_v = 360
                shift_start = (max(v1)+min(v2))/2
                new_values = map(lambda x: x+shift_v if x<shift_start else x, values)
                values = sorted(new_values)
        # limits for 100%, 99.99% 99.9%, 99.5%, 99%, 95%
        for pr in [100.0, 99.99, 99.9, 99.5, 99.0, 95.0]:
            nn = int(math.ceil(n*pr/100.0))
            assert nn>0
            if key in ['nn_ang_norm']:
                i1 = 0
            elif key in ['n2_z']:
                if min(values)<-0.5:
                    i1 = 0
                elif max(values)>0.5:
                    i1 = n-nn
            else:
                i1 = (n-nn)//2
            i2 = i1+nn
            assert len(values[i1:i2])==nn
            lost = n-nn
            limit_key = "%s/%s/%s/%.2f" % (options.desc, options.n_type, key, pr)
            min_v = values[i1]
            max_v = values[i2-1]
            if key in ['rot_ang'] and max_v>270.0:
                l = {'min1': min_v, 'max1': 270.0, 'min2': -90.0, 'max2': max_v-shift_v}
            else:
                l = {'min': min_v, 'max': max_v}
            limit_str = "limits['%(limit_key)s'] = %(l)s # lost: %(lost)d" % locals() 
            all_limits_str += limit_str + "\n"
            print >> sys.stderr, "key",key,"pr",pr,"min",values[i1],"max",values[i2-1],"shift",shift_v,"lost",lost,"limit",l
    out_fn = os.path.join(options.output_dir,"%s_%s.txt" % (options.desc,options.n_type))
    ensure_dir(out_fn)
    write_file(out_fn, all_limits_str)


def main():
    (_parser, options, _args) = parse_args()
    
    gen_limits(options)

if __name__=="__main__":
    main()