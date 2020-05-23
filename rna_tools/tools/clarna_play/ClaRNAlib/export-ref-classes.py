#!/usr/bin/env python
#
# TODO: dodać konwersję nazw ze znormalizowanych na oryginalne!!!
#
import sys
import os
import re
import math
import gzip
import zipfile
from optparse import OptionParser
import xlwt

from utils import load_json, DoubletDescTool, N_TYPES

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""export reference classes to XLS file""")
    parser.add_option("-i", dest="input",
                  help="read groups from file", metavar="FILE",
                  default=FileNamesObject.groups_fn(setname="training",reduced=True))
    parser.add_option("--output-xls", dest="output_xls",
                  help="save result to output XLS", metavar="FILE")
    parser.add_option("--output-zip", dest="output_zip",
                  help="save result to output ZIP", metavar="FILE")
    parser.add_option("--xls-rows-limit", dest="xls_rows_limit",
                  help="limit on number of rows in XLS files", metavar="N", default=10000)
    parser.add_option("--txt-rows-limit", dest="txt_rows_limit",
                  help="limit on number of rows in TXT files", metavar="N", default=1000000)

    (options, args)  = parser.parse_args()
    
    options.xls_rows_limit = int(options.xls_rows_limit)
    options.txt_rows_limit = int(options.txt_rows_limit)
    return (parser, options, args)

def export_classes_to_xls(options):
    data = load_json(options.input)
    
    wb = xlwt.Workbook()
    
    n_types = sorted(N_TYPES,key=DoubletDescTool.n_type_sort_key)
    
    for sc in ['bp','stacking','base-phosphate','base-ribose','other','other2','other3']:
        descs = [x.split("/")[2] for x in data.keys() if re.match("^classifier/%s/"%sc,x)]
        descs = sorted(list(set(descs)),key=DoubletDescTool.desc_sort_key)

        for desc in descs:
            title = 'Reference groups for: %s/%s'%(sc,desc)
            print "adding %s" % title
            ws = wb.add_sheet(sc+"_"+desc)
            ws.write(0,0,title)

            ws.write(2,0,"N. Type")
            ws.write(3,0,"Count")
            for i,n_type in enumerate(n_types):
                ws.write(2,1+i,n_type)
                k = "classifier/%s/%s/%s" % (sc,desc,n_type)
                elems = data.get(k)
                elems.sort()
                if len(elems)>options.xls_rows_limit:
                    print "reducing size of %s from %d to %d" % (k,len(elems),options.xls_rows_limit)
                    elems = elems[0:options.xls_rows_limit]
                ws.write(3,1+i,len(elems))
                for j,d_id in enumerate(sorted(elems)):
                    ws.write(4+j,1+i,d_id)
    wb.save(options.output_xls)

def export_classes_to_zip(options):
    data = load_json(options.input)

    n_types = sorted(N_TYPES,key=DoubletDescTool.n_type_sort_key)
    
    with zipfile.ZipFile(options.output_zip, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for sc in ['bp','stacking','base-phosphate','base-ribose','other','other2']:
            descs = [x.split("/")[2] for x in data.keys() if re.match("^classifier/%s/"%sc,x)]
            descs = sorted(list(set(descs)),key=DoubletDescTool.desc_sort_key)
    
            for desc in descs:
                title = 'Reference groups for: %s/%s'%(sc,desc)
                print "adding %s" % title

                safe_desc = desc.replace("<","lt").replace(">","gt")
                for i,n_type in enumerate(n_types):
                    k = "classifier/%s/%s/%s" % (sc,desc,n_type)
                    elems = data.get(k)
                    if elems is not None:
                        elems.sort()
                        if len(elems)>options.txt_rows_limit:
                            print "reducing size of %s from %d to %d" % (k,len(elems),options.txt_rows_limit)
                            elems = elems[0::options.txt_rows_limit]
                        z.writestr("ref/%(sc)s/%(safe_desc)s/%(sc)s_%(safe_desc)s_%(n_type)s.txt"%locals(), "\n".join(elems))

def main():
    (parser, options, args) = parse_args()
    if options.output_xls:
        export_classes_to_xls(options)
    elif options.output_zip:
        export_classes_to_zip(options)
    else:
        raise Exception("Select output mode")
    
if __name__ == '__main__':
    main()