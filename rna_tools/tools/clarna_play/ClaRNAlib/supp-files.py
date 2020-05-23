#!/usr/bin/env python
import re
import subprocess
import urllib
from optparse import OptionParser
from utils import *
import xlwt

import numpy as np

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""evaluation graphs""")
    parser.add_option("--output-dir", dest="output_dir",
                  help="output PNG", metavar="DIR", default="doc/nar-supp")
    parser.add_option("--force", dest="force", action='store_true', 
                  help="force creation of files",default=False)

    (options, args)  = parser.parse_args()
    return (parser, options, args)

def cl_name(x):
    if x=="CL":
        return "ClaRNA"
    elif x=='CL_ignore_fuzzy':
        return "ClaRNA (with ignored fuzzy doublets)"
    elif x=="RV":
        return "RNA-View"
    elif x=="MC":
        return "MC-Annotate"
    elif x=="MO":
        return "ModeRNA"
    elif x=="FR":
        return "FR3D"
    else:
        return "UNKNOWN"

class ExcelWriter():

    def __init__(self):
        font_bold = xlwt.Font()
        font_bold.name = 'Arial'
        font_bold.bold = True

        style_norm = xlwt.XFStyle()
        style_bf = xlwt.XFStyle()
        style_bf.font = font_bold
        
        self.styles = {}
        self.styles['normal'] = style_norm
        self.styles['normal_bf'] = style_bf
        self.styles['bf'] = style_bf
        self.styles['hdr'] = style_bf
        self.styles['num'] = xlwt.easyxf('',num_format_str='0.0000')
        self.styles['num_bf'] = xlwt.easyxf('font: bold on',num_format_str='0.0000')
        
        self.wb = xlwt.Workbook()
        self.ws = None
        
    def add_sheet(self,sheet_name):
        self.ws = self.wb.add_sheet(sheet_name)
    
    def write(self,row,col,cellText,style='normal'):
        assert self.ws is not None
        self.ws.write(row,col,cellText,self.styles.get(style))

    def save(self,fn):
        self.wb.save(fn)

##############################################
    
def supp_dataset_pdb_list(options):
    status_msg("generating supp-dataset-pdb-list")

    fn = os.path.join(options.output_dir,"supp-file-s1-list-of-pdb-files.xls")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    
    fields = ["pdb_id","title","release_date","resolution","url"]
    header = ["PDB ID","Title","Release date","Resolution","URL"]
    
    e = ExcelWriter()
    for sheet_name,data_fn in [
            ('Training Set',"data-training.txt"),
            ('Testing Set',"data-bench.txt")
        ]:
        e.add_sheet(sheet_name)
        e.write(0,0,sheet_name,style='bf')
        for i,h in enumerate(header):
            e.write(2,i,h,style='bf')
        for i,pdb_id in enumerate(read_file(data_fn).strip().split("\n")):
            print " - adding %s" % pdb_id
            row = get_pdb_info(pdb_id)
            print " - row: %s" % row
            if not re.match("^zz",pdb_id):
                row["url"]="http://www.pdb.org/pdb/explore/explore.do?structureId="+pdb_id
            for j,field in enumerate(fields):
                e.write(3+i,j,row.get(field,""))
        e.ws.col(1).width = 256*50
    e.save(fn)


def supp_fig_clarna_workflow(options):
    status_msg("generating supp-fig-clarna-workflow")

    fn = os.path.join(options.output_dir,"supp-figure-s1-clarna-workflow.pdf")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    r1 = os.system("cd doc && make nar-supp-fig-workflow.pdf");
    assert r1==0
    r2 = os.system("pdfcrop doc/nar-supp-fig-workflow.pdf '%s'" % fn)
    assert r2==0

def supp_fig_bph(options):
    status_msg("generating supp-fig-bph")

    fn = os.path.join(options.output_dir,"supp-figure-s2-bph.png")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    environ = dict(os.environ)
    environ['SHELL'] = '/bin/bash'
    r1 = subprocess.call("fab-2.7 nar_images_bph",shell=True,env=environ)
    assert r1==0
    r2 = os.system("cp 'doc/supp-figure-2-bph.png' '%s'" % fn)
    assert r2==0


def supp_fig_clarna_time(options):
    status_msg("generating supp-fig-clarna-time")

    fn = os.path.join(options.output_dir,"supp-figure-s4-clarna-time-benchmarks.png")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    environ = dict(os.environ)
    environ['SHELL'] = '/bin/bash'
    r1 = subprocess.call("fab-2.7 nar_images_benchmarks",shell=True,env=environ)
    assert r1==0
    r2 = os.system("cp 'doc/nar-fig-clarna-time2.png' '%s'" % fn)
    assert r2==0

def _supp_table_classifiers_eval_params(e,options):
    data = load_json(FileNamesObject.eval_fn(setname="bench",t="eval3"))

    e.add_sheet("Classifiers evaluation")
    e.write(0,0,"Classifiers evaluation","bf")
    
    row = 3

    fields = ["tp","tn","fp","fn","tpr","fpr","acc","spc","ppv","npv","fdr","mcc","f1"]

    col = len(fields)+2
    e.write(row,col,"Classifier codes","bf")
    # for i,prg in enumerate(["CL","CL_ignore_fuzzy","FR","MC","MO","RV"]):
    for i,prg in enumerate(["CL","FR","MC","MO","RV"]):
        e.write(row+1+i,col,prg)
        e.write(row+1+i,col+1,cl_name(prg))
    
    for sc in ['bp-classic','bp-non-classic','stacking','base-phosphate','base-ribose']:
        e.write(row,0,sc,"bf")
        row += 1
        for j,f in enumerate(fields):
            e.write(row,1+j,f,'bf')
        row += 1
        if sc in ['bp-classic','bp-non-classic']:
            programs = ['CL','CL_ignore_fuzzy','FR','MC','RV']
        elif sc in ['stacking']:
            programs = ['CL','CL_ignore_fuzzy','FR','MC','MO']
        else:
            programs = ['CL','CL_ignore_fuzzy']
        programs.remove('CL_ignore_fuzzy')

        res = {}
        for i,prg in enumerate(programs):
            r = {}
            for k in ['tp','fp','tn','fn']:
                r[k] = data.get("evaluation/%(prg)s/%(k)s/%(sc)s" % locals(),0)
            r = confusion_matrix_params(r)
            res[prg] = r
        best = {}
        for k in res.values()[0].keys():
            if k in ['fp','fn','fpr','fdr']:
                best[k] = min([res[prg][k] for prg in programs])
            else:
                best[k] = max([res[prg][k] for prg in programs])

        for i,prg in enumerate(programs):
            r = res[prg]
            e.write(row,0,prg,"bf")
            for j,f in enumerate(fields):
                v = r[f]
                if f in ['tp','fp','tn','fn']:
                    s = 'normal'
                else:
                    s = 'num'
                if len(programs)>1 and v==best[f]:
                    s += '_bf'
                e.write(row,1+j,v,s)
            row += 1
        
        row += 1

def _supp_table_classifier_5_cross_validation(e,options):
    return _supp_table_classifier_cross_validation("Results from 5-cross-validation",range(1,6),e,options)

def _supp_table_classifier_13_cross_validation(e,options):
    return _supp_table_classifier_cross_validation("Results from 13-cross-val.",[chr(ord('A')+i) for i in range(13)],e,options)

def _supp_table_classifier_cross_validation(title,nums,e,options):
    import tarfile
    tmpdir = "/tmp/cross_results"
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)
    
    results = {}
    for num in nums:
        tar = tarfile.open("gc-data/cross-%s.tgz"%num,mode="r:gz")
        stats_fn = "eval/"+os.path.basename(FileNamesObject.eval_fn(setname="training",t="eval3"))
        tar.extract(stats_fn,tmpdir)
        data = load_json(os.path.join(tmpdir,stats_fn))
        for sc in ['bp-classic','bp-non-classic','stacking','base-phosphate','base-ribose']:
            r = {}
            for k in ['tp','fp','tn','fn']:
                r[k] = data.get("evaluation/CL/%(k)s/%(sc)s" % locals(),0)
            r = confusion_matrix_params(r)
            if not results.has_key(sc):
                results[sc] = {}
            for k,v in r.items():
                if not results[sc].has_key(k):
                    results[sc][k] = []
                results[sc][k].append(v)

    e.add_sheet(title)
    e.write(0,0,title,"hdr")
    e.write(1,0,"average results over %d cases"%len(nums),"normal")
    
    row = 3
    
    fields = ["tp","tn","fp","fn","tpr","fpr","acc","spc","ppv","npv","fdr","mcc","f1"]

    for sc in ['bp-classic','bp-non-classic','stacking','base-phosphate','base-ribose']:
        e.write(row,0,sc,"bf")
        row += 1
        for j,f in enumerate(fields):
            e.write(row,1+j,f,'bf')
        row += 1
        for op in ["avg","std-dev","min","max"]:
            r = results[sc]
            e.write(row,0,"CL (%s)"%op,"bf")
            for j,f in enumerate(fields):
                print sc, f, np.mean(r[f]), min(r[f]), max(r[f])
                if op=='avg':
                    v = np.mean(r[f])
                elif op=='std-dev':
                    v = np.std(r[f])
                elif op=='min':
                    v = min(r[f])
                elif op=='max':
                    v = max(r[f])
                else:
                    raise Exception("unknown op: %s"%op)

                if f in ['tp','fp','tn','fn']:
                    v = int(round(v))
                    s = 'normal'
                else:
                    s = 'num'
                e.write(row,1+j,v,s)
            row += 1
        row += 1

def _supp_table_classifier_similarity(e,options):
    e.add_sheet("Similarity between classifiers")
    e.write(0,0,"Similarity between classifiers","hdr")
    e.write(2,0,"The similarity score is obtained using formula:")
    e.write(3,0,"number of doublets recognized by both classifiers with the same interaction type divided by number of doublets recognized by any of the classifiers")
    
    
    data = load_json(FileNamesObject.eval_fn(setname="bench",t="corr"))
    row = 5

    e.write(row,6,"Classifier codes","bf")
    for i,prg in enumerate(["CL","FR","MC","MO","RV"]):
        e.write(row+1+i,6,prg)
        e.write(row+1+i,7,cl_name(prg))

    for sc in ['bp-classic','bp-non-classic','stacking','base-phosphate','base-ribose']:
        e.write(row,0,sc,"bf")
        row += 1
        if sc in ['bp-classic','bp-non-classic']:
            programs = ['CL','FR','MC','RV']
        elif sc in ['stacking']:
            programs = ['CL','FR','MC','MO']
        else:
            programs = ['CL','FR']
        for i,prg in enumerate(programs):
            e.write(row,1+i,prg,"bf")
            e.write(row+1+i,0,prg,"bf")
        for i,prg1 in enumerate(programs):
            for j,prg2 in enumerate(programs):
                p1 = min(prg1,prg2)
                p2 = max(prg1,prg2)
                all = data.get('correlation/%(sc)s/%(p1)s_vs_%(p2)s/all'%locals(),0)
                common = data.get('correlation/%(sc)s/%(p1)s_vs_%(p2)s/common'%locals(),0)
                v = float(common)/max(all,1)
                if prg1==prg2:
                    v = 1.000
                e.write(row+1+i,1+j,v,'num')
        
        row += len(programs)+3

def supp_table_classifiers_evaluation_params(options):
    status_msg("generating supp-table-classifiers-evaluation")
    
    fn = os.path.join(options.output_dir,"supp-table-s1-classifiers-evaluation.xls")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return
        
    e = ExcelWriter()

    # wyniki klasyfikatorow w tym clarny na testingset
    _supp_table_classifiers_eval_params(e,options)
    # wyniki similarity score dla wszystkich klasyfikatorow dla testing set
    _supp_table_classifier_similarity(e,options)
    # wyniki (average) dla 5-cross validation on training set
    _supp_table_classifier_5_cross_validation(e,options)
    _supp_table_classifier_13_cross_validation(e,options)

    e.save(fn)


def supp_table_results_diff(options):
    status_msg("generating supp-table-clarna-results")
    
    fn = os.path.join(options.output_dir,"supp-table-web-only-results-differences.xls")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return
        
    data_t = load_json(FileNamesObject.eval_fn(setname="training",t="eval3"))
    data_v = load_json(FileNamesObject.eval_fn(setname="bench",t="eval3"))

    e = ExcelWriter()
    
    for suffix,sheet_name,sc_list in [
        ('new-base-ribose','New Base Ribose',['base-ribose']),
        ('cis-vs-trans','Differences on Cis and Trans',['bp']),
        ('wh_cis-vs-sh_cis','WH vs SH_cis',['bp']),
        #('consistent-with-single-cl','Consistent with single classifier',['bp','stacking','base-phosphate','base-ribose']),
        #('not-recognized-by-others-cl','Not detected by other classifiers',['bp','stacking','base-phosphate','base-ribose']),
        #('others','Possibly wrong detection',['bp','stacking','base-phosphate','base-ribose'])
        ]:
        e.add_sheet(sheet_name.replace(" ","_")[0:16])
        e.write(0,0,sheet_name,"hdr")

        row = 2
        for sc in sc_list:
            if sc in ['bp']:
                PRGS = ["FR","MC","RV"]
            elif sc in ['stacking']:
                PRGS = ["FR","MC","MO"]
            elif sc in ['base-phosphate','base-ribose']:
                PRGS = ["FR"]
            else:
                raise Exception("bad sc: %s"%sc)

            fields = ['set_name','pdb_id','res1','res2','exp_res']+['res_'+x for x in PRGS]+['cl_res','extra']
            header = ['Set name','PDB','Residue-A','Residue-B','Expected result']+[x+" result" for x in PRGS]+['Clarna result','Extra information']
            e.write(row,0,sc,'bf')
            row += 1
            for i,h in enumerate(header):
                e.write(row,i,h,'bf')
            row += 1
            
            count = 0
            data_key = 'evaluation-details/CL/fp-%s/%s'%(suffix,sc)
            if sc=='stacking' and suffix=='consistent-with-single-cl':
                e.write(row,0,"%s examples in training set" % len(data_t.get(data_key)))
                row += 1
                e.write(row,0,"%s examples in testing set" % len(data_v.get(data_key)))
                row += 1
                e.write(row,0,"details available on ClaRNA web-page")
                row += 1
            else:
                for set_name,data in ('Training',data_t),('Testing',data_v):

                    for r in sorted(data.get(data_key,[]),key=lambda x:x['d_id']):
                        if row>60000:
                            continue
                        print r
                        r['set_name'] = set_name
                        r['pdb_id'],r['res1'],r['res2'] = r['d_id'].split(":")
                        r['pdb_id'] = r['pdb_id'].lower()
                        for p,o in r['other_results']:
                            r['res_'+p] = o
                        for i,f in enumerate(fields):
                            e.write(row,i,r.get(f,''),'normal')
                        row += 1
                        count += 1
                if count==0:
                    e.write(row,0,"no examples")
                    row += 1
            row += 1
    e.save(fn)

######################################
## Web-Only files
######################################

def supp_table_clarna_results(options):
    status_msg("generating supp-table-clarna-results")
    
    url = "http://rnacontacts-vm/alg/new-classifier-evaluation/xls/bench"
    fn = os.path.join(options.output_dir,"supp-table-web-only-clarna-results.xls")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return
    r = os.system("wget '%(url)s' -O '%(fn)s'" % locals())
    assert r==0

def supp_fig_ref_doublets(options):
    status_msg("generating supp-fig-ref-doublets")

    fn = os.path.join(options.output_dir,"supp-figure-web-only-reference-doublets.pdf")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    status_msg(" *** NOT IMPLEMENTED YET! *** ")

def supp_table_descriptions_dict(options):
    status_msg("generating supp-table-interaction-types-dictionary")

    fn = os.path.join(options.output_dir,"supp-table-web-only-interaction-types-dictionary.xls")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    d = load_json("descriptions-dict.json")

    e = ExcelWriter()
    for sheet_name,prg in [
            ('RNA-View','RV'),
            ('MC-Annotate','MC'),
            ('FR3D','FR'),
        ]:
        e.add_sheet(sheet_name)
        
        e.write(0,0,sheet_name+" - dictionary","bf")

        header = ["Name used by %s"%sheet_name,"Name used by ClaRNA"]
        
        r_num = 1
        for sc in ['bp','stacking','base-ribose','base-phosphate']:
            elems = [(k,v) for k,v in sorted(d[prg].items()) if DoubletDescTool.get_desc_category(v)==sc]
            if len(elems)>0:
                r_num += 2
                e.write(r_num,0,"Interaction class: "+sc,"bf")
                r_num += 1
                for i,h in enumerate(header):
                    e.write(r_num,i,h,"bf")

                for k,v in elems:
                    k = re.sub("^stacking_","",k)
                    r_num += 1
                    e.write(r_num,0,k)
                    e.write(r_num,1,v)
        if prg in ['MC']:
            e.ws.col(0).width = 256*50
        else:
            e.ws.col(0).width = 256*20
        e.ws.col(1).width = 256*20

    e.save(fn)

def supp_dataset_reference_set(options):
    status_msg("generating supp-dataset-ref-set")

    fn = os.path.join(options.output_dir,"supp-dataset-web-only-reference-set.zip")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    status_msg(" *** NOT IMPLEMENTED YET! *** ")


def _removed_supp_table_changes_in_ref_set(options):
    # usuniete zgodnie z sugestia JMB (2013--08-09)
    status_msg("generating supp-table-changes-in-ref")
    
    fn = os.path.join(options.output_dir,"supp-table-web-only-list-of-changes-in-training-set.xls")
    if os.path.isfile(fn) and not options.force:
        status_msg(" - skipping generation of %s, file already exists" % fn)
        return

    not_ref = {}
    for x in find_files_in_dir("expert/not-ref"):
        if re.match("^.*.json$",x):
            tmp = x.split(".")[0].split("_")
            assert len(tmp)>=4
            sc = tmp[0]
            n_type = tmp[1]
            desc = tmp[2]+"_"+tmp[3]
            if len(tmp)>4:
                reason = " ".join(tmp[4:])
            else:
                reason = ""
            
            for d_id in load_json(x):
                if not_ref.has_key(d_id):
                    print "WARNING: %s has been defined twice! %s,%s" % (d_id,not_ref[d_id],x)
                not_ref[d_id] = {"sc":sc,"n_type":n_type,"exp_desc":desc,"reason":reason}

    e = ExcelWriter()
    for sheet_name,data_fn in [
            ('Changes in Training Set',"data-training.txt"),
        ]:
        e.add_sheet(sheet_name)
        
        pdb_ids = set(read_file(data_fn).strip().split("\n"))
        
        e.write(0,0,sheet_name,"bf")

        fields = ["pdb_id","residue_a","residue_b","n_type","exp_desc","reason"]
        header = ["PDB ID","Residue A","Residue B","N Type","Expected result","Reason of exclusion from the reference set"]
        for i,h in enumerate(header):
            e.write(2,i,h,"bf")
        
        r_num = 2
        for d_id in sorted(not_ref.keys()):
            pdb_id,r1,r2 = d_id.split(":")
            pdb_id = pdb_id.lower()
            row = not_ref[d_id]
            row['pdb_id'] = pdb_id
            row['residue_a'] = r1
            row['residue_b'] = r2
            if pdb_id in pdb_ids:
                r_num += 1
                for j,f in enumerate(fields):
                    e.write(r_num,j,row.get(f))

    e.save(fn)

######################################

def add_labels():
    for fn,label,pointsize,dx,dy in [
        ("supp-figure-2-bph.png","Supplementary Figure 2. Example of 3 close doublets from merged BPh classes (W_3/4/5BPh)",32,100,10),
        ("supp-figure-3-merged-bph-classes.png","Supplementary Figure 3. Clusters of doublets from merged BPh classes (W_3/4/5BPh and H_7/8/9BPh)",50,200,10),
        ("supp-figure-4-clarna-time-benchmarks.png","Supplementary Figure 4. Time efficiency benchmarks of ClaRNA",50,10,10),
    ]:
        full_fn = os.path.join("doc","nar-supp",fn)
        out_fn = "/tmp/"+fn
        cmd = "montage -background white -font Arial -pointsize %(pointsize)d -label '\n%(label)s' '%(full_fn)s' -geometry '+%(dx)d+%(dy)d' '%(out_fn)s'" % locals()
        print cmd
        res = os.system(cmd)
        assert res==0

def main():
    (parser,options,_args) = parse_args()
    options.force = True
    supp_table_classifiers_evaluation_params(options)
    # supp_table_results_diff(options)
    #or
    # add_labels()
    return
    
    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)
    
    supp_dataset_pdb_list(options)
    supp_dataset_reference_set(options)
    supp_fig_clarna_time(options)
    supp_fig_clarna_workflow(options)
    supp_fig_ref_doublets(options)
    supp_fig_bph(options)
    # TODO zautomatyzowac 
    # fab tw_20130619_nar_supp_fig_merged_bph_classes
    supp_table_classifiers_evaluation_params(options)
    supp_table_results_diff(options)

    # web-only
    # supp_table_changes_in_ref_set(options)
    supp_table_descriptions_dict(options)
    supp_table_clarna_results(options)

if __name__=="__main__":
    main()
