#!/usr/bin/env python
import re
from optparse import OptionParser
from utils import *
from distances import bp_distance1, bp_distance2, bp_distance3, expected_strand_orient

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""evaluation graphs""")
    parser.add_option("-i", "--input", dest="input",
                  help="input JSON", metavar="JSON", default=FileNamesObject.eval_fn(setname="bench",t="eval3"))
    parser.add_option("-o", "--output-png", dest="output_png",
                  help="output PNG", metavar="PNG")
    parser.add_option("--output-pdf", dest="output_pdf",
                  help="output PDF", metavar="PDF")
    (options, args)  = parser.parse_args()
    return (parser, options, args)

def main():
    (parser,options,_args) = parse_args()
    
    data = load_json(options.input)
    
    DRAW_TWO_PLOTS = False
    
    sub_categories = ['bp-classic','bp-non-classic','stacking','base-phosphate','base-ribose']
    sub_categories_labels = ['classical pairs','non-classical pairs','stackings','base-phosphate','base-ribose']
    
    eps = 0.00001

    matplotlib.rcParams.update({'font.size': 22})
    
    if DRAW_TWO_PLOTS:
        f, (ax,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(20,7))
    
        ax.margins(0.1)
        ax2.margins(0.1)
    else:
        f = plt.figure(figsize=(10,9))
        ax = f.add_subplot(111)
        ax.margins(0.1)
    marks = ["o","o","s","v","^","*"]
    colors = ["red","green","blue","orange","violet","cyan"]
    for sc_num,(sc,sc_label) in enumerate(zip(sub_categories,sub_categories_labels)):
        keys = sorted([x.split("/")[-1] for x in data.keys() if re.match('^evaluation/CL/tp/'+sc+"/[0-9.]+$",x)],key=lambda x: float(x))
        print keys
        if len(keys)==0:
            continue
        roc_data = []
        mark_data = []
        roc_data2 = []
        for k in keys:
            if sc=='base-ribose' and float(k)<0.1:
                continue
            row = {}
            for i in ['tp','fp','fn','tn']:
                row[i] = data.get("evaluation/CL/"+i+"/"+sc+'/'+k,0)
            row = confusion_matrix_params(row)
            row["sc"] = sc
            row["k"] = float(k)
            print "%(sc)-20s k=%(k)4.1f tp=%(tp)7d tn=%(tn)7d fp=%(fp)7d fn=%(fn)7d tpr=%(tpr).3f fpr=%(fpr).3f acc=%(acc).3f fdr=%(fdr).3f" % row
            # print sc,k,"tpr",row['tpr'],"fpr",row['fpr'],"acc",row["acc"],"fdr",row["fdr"]
            roc_data.append((row['fpr'],min(1.0-eps,row['tpr'])))
            if float(k)>=0.0 or True:
                roc_data2.append((row['fpr'],row['tpr']))
            if float(k)==0.5:
                mark_data.append((row['fpr'],min(1.0-eps,row['tpr'])))

        c = colors[sc_num%len(colors)]
        ax.plot([x[0] for x in roc_data], [1.0-x[1] for x in roc_data], marks[sc_num%len(marks)]+'-', label=sc_label, ms=10, color=c, linewidth=3.0)
        if DRAW_TWO_PLOTS:
            ax2.plot([x[0] for x in roc_data2], [x[1] for x in roc_data2], marks[sc_num%len(marks)]+'-', label=sc_label, color=c)

        ax.plot([x[0] for x in mark_data], [1.0-x[1] for x in mark_data], marks[sc_num%len(marks)], ms=25, color=c, linewidth=3.0)
        if DRAW_TWO_PLOTS:
            ax2.plot([x[0] for x in mark_data], [x[1] for x in mark_data], marks[sc_num%len(marks)], ms=15, color=c)

    ax.plot(np.arange(0.0,0.2,0.01),[1.0-x for x in np.arange(0.0,0.2,0.01)],':', color="gray",linewidth=3.0)
    ax.text( 0.15, 0.7, 'TPR=FPR', color="gray")

    # plt.title(sc)
    if DRAW_TWO_PLOTS:
        ax.set_title("(a)")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.tick_params(axis='y', pad=8)

    labelx = -0.15
    ax.yaxis.set_label_coords(labelx, 0.5)
    ax.set_xlim(-0.01,0.2)
    ax.set_ylim(0.0,1.5)

    if DRAW_TWO_PLOTS:
        ax2.set_title("(b)")
        ax2.set_ylabel("True positive rate")
        ax2.set_xlabel("False positive rate")
        ax2.yaxis.set_label_coords(-0.1, 0.5)
        ax2.set_ylim(0.95,1.005)

    ax.set_yscale('log',basey=10)
    ax.invert_yaxis()
    print ax.get_ylim()
    ax.set_ylim(1.5,ax.get_ylim()[1]*0.5)
    # plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    
    def one_minus(tick_val, tick_pos):
        if tick_val>=0.01:
            return "%.2f" % (1-tick_val)
        if tick_val>=0.001:
            return "%.3f" % (1-tick_val)
        v = math.log(tick_val,10)
        # print "ONE_MINUS", v, int(round(v))
        return '$1-10^{%d}$' % int(round(v))

    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(one_minus)) 

    plt.subplots_adjust(left=0.17, right=0.96, top=0.99, bottom=0.1)

    if DRAW_TWO_PLOTS:
        plt.legend(loc='lower right')
    else:
        plt.legend(loc=(0.51,0.15),prop={'size':20})
    if options.output_png:
        plt.savefig(options.output_png)
    elif options.output_pdf:
        plt.savefig(options.output_pdf,dpi=300)
    else:
        plt.show()

if __name__=="__main__":
    main()
