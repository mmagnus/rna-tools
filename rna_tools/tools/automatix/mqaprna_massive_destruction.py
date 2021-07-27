#!/usr/bin/python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
1 - method
2 - cpu
3 - "_X.csv' (without csv

   print "1 - method, 2 - cpu, 3 - '_X.csv' (without csv) 4 - run 5 - method_code, e.g. Fhi"
"""
from __future__ import print_function

import argparse
import sys
import os
import glob


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("method")
    parser.add_argument("cpu", type=int)
    parser.add_argument("--run", action="store_true") # nargs='+')
    parser.add_argument("path")
    parser.add_argument("--subfolder", default="struc", help="struc in this case /home/mqapRNA/mqaprna_datasets/rasp/2nueC/struc/2nueC_M10.pdb")
    parser.add_argument("dirs", nargs='+')# help="", default="") # 
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    missing = ''
    
    files = os.listdir(PATH)
    method = args.method
    cpu = int(args.cpu)
    method_code = ''#sys.argv[5]

    for f in args.dirs:
        if os.path.isdir(f):
            print('#' * 80)
            print(f)
            if missing:
                if f not in  missing.split(','):
                    continue

            #cmd = 'cp -v ' + f + '/' + f + '.pdb ' + f + '/struc/'
            #print cmd
            #cmd = "mqaprna.py -n 1a9nR/1a9nR.pdb 1a9nR/struc/ 1a9nR"
            #cmd = "mqaprna.py -m 6 -t SimRNA_1 " + f + '/struc ' + f + '_simrna_1_'

            #cmd = "mqaprna.py -m " + str(cpu) + " -g /home/magnus/mqaprna_decoys/rasp/" + f + '_' + method + "-v0.3-48-g80cbfc0-dirty.csv -t " + method + "  " + PATH + f + '/struc ' + PATH + '_' + f + '_' + method 

            pattern = path + '/*' + f + '*' + method + '**csv'
            print(pattern)
            csvfn_to_ignore = glob.glob(pattern)
            print(csvfn_to_ignore)
            if len(csvfn_to_ignore) > 0:
                g = " -g " + csvfn_to_ignore[0]
            else:
                print('nothing to ignore...')
                g = ''
            #g = '' # g hack!

            if cpu > 0:
                multi_cpu = ' -m ' + str(cpu)
                cpu_c = cpu
            else:
                multi_cpu = ''
                cpu_c = '1'

            cmd = "mqaprna.py " + multi_cpu + " "  + g + " -t " + method + "  " + PATH + f + os.sep + args.subfolder + '/*.pdb -o ' + PATH + '' + f + '_' + method
            print(cmd)

            #cmd2 = 'echo "' + cmd + '" | qsub -cwd -V -pe mpi ' + str(cpu_c) + '  -N _' + f + method_code #+ f + '_'+method #  -l h_vmem=300M -l mem_free=500M
            # print cmd2
            #RASP,SimRNA,AnalyzeGeometry,FARNA,QRNA,NAST_pyro,radiu
            #                s_of_gyration,SSAgreement,ClashScore,RNAkb
            if args.run:
                os.system(cmd)

