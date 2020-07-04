#!/usr/bin/python
# -*- coding: utf-8 -*-
"""rna_clanstix - a tool for visualizing RNA 3D structures based on pairwise structural similarity with Clans.

We hacked Clans thus instead of BLAST-based distances between sequences, you can analyze distances between structures described as p-values of rmsd (based on the method from the Dokholyan lab.)

Quickref::

      rna_clanstix.py --groups-auto 10 --color-by-homolog --shape-by-source  thf_ref_mapping_pk_refX.txt input2.clans

Running Clans:
To run CLANS you need to have Java 1.4 or better installed (java can be downloaded HERE). For full functionality you will also need the NCBI BLAST,PSI-BLAST and formatdb executables (NCBI). For command line parameters and basic help please refer to the README file.
(source: http://www.eb.tuebingen.mpg.de/research/departments/protein-evolution/software/clans.html)

.. image:: ../../rna_tools/tools/clanstix/doc/yndSrLTb7l.gif

The RMSDs between structures are converted into p-values based on the method from the Dokholyan lab or some hacky way developed by mmagnus .

Color groups
---------------------------------------
You can color your groups:

.. image:: ../../rna_tools/tools/clanstix/doc/rna_clanstix.png

To get colors, run a cmd like this::

   rna_clastix.py rnapz17_matrix_farfar_HelSeedCst.txt --groups 20:seq1+20+20+20+20+20+20:seq10

where with the ``+`` sign you separate groups. Each group has to have a number of structures. Optionally it can have a name, e.g., ``20:seq1``, use ``:`` as a separator. If a provided name is ``native`` then this group will be shown as starts.

Get inspiration for more colors (http://www.rapidtables.com/web/color/RGB_Color.htm)

How to use ClanstixRNA?
----------------------------------------

1. Get a matrix of distances, save it as e.g. matrix.txt (see Comment below)
2. run ClanstixRNA on this matrix to get an input file to Clans (e.g. clans_rna.txt)::

     rna_clanstix.py test_data/matrix.txt # clans.input will be created by default

3. open CLANS and click File -> Load run and load clans_run.txt
4. You're done! :-)

Comment: To get this matrix you can use for example another tool from the rna-pdb-tools packages::

     rna_calc_rmsd_all_vs_all.py -i rp18 -o rp18_rmsd.csv
     rna_clastix.py --groups 1:native+5:3dRNA+
           5:Chen+3:Dokh+5:Feng+5:LeeASModel+
           5:Lee+5:RNAComposer+10:RW3D+5:Rhiju+
           1:YagoubAli+3:SimRNA  rp18_rmsd.csv clans.in

     rna_clastix.py --groups 100+100+100+100+100+100+100+100+100+100+1:native  rp18_rmsd.csv

where ``rp18`` is a folder with structure and ``rp18_rmsd.csv`` is a matrix of all-vs-all rmsds.

.. image:: ../../rna_tools/tools/clanstix/doc/rp18_clanstix.png

Hajdin, C. E., Ding, F., Dokholyan, N. V, & Weeks, K. M. (2010). On the significance of an RNA tertiary structure prediction. RNA (New York, N.Y.), 16(7), 1340–9. doi:10.1261/rna.1837410

An output of this tool can be viewed using CLANS.

Frickey, T., & Lupas, A. (2004). CLANS: a Java application for visualizing protein families based on pairwise similarity. Bioinformatics (Oxford, England), 20(18), 3702–4. doi:10.1093/bioinformatics/bth444
"""
from __future__ import print_function
import argparse
import rna_tools.tools.rna_prediction_significance.rna_prediction_significance as pv
import numpy as np
import math
import logging
import time

logging.basicConfig(level=logging.INFO,
                format='%(message)s',
                datefmt='%m-%d %H:%M',
                filename='rna_clanstix.log',
                filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

log = logging.getLogger()
log.setLevel(logging.INFO)


class RNAStructClans:
    """Clans run.

    Usage::

        >>> f = open('matrix.txt')
        >>> ids = f.readline().replace('#','').split()
        >>> c = RNAStructClans(n=len(ids)) # 200?
        >>> c.add_ids(ids)
        >>> c.dist_from_matrix(f)
        >>> print(c.txt)
    """

    def __init__(self, n=10, dotsize=10):
        self.n = n
        self.comment = ''
        #cluster2d=false
        self.txt = """sequences=%i
<param>
maxmove=0.1
pval=%s
usescval=false
complexatt=true
cooling=1.0
currcool=1.0
attfactor=10.0
attvalpow=1
repfactor=10.0
repvalpow=1
dampening=1.0
minattract=1.0
cluster2d=true
blastpath=''
formatdbpath=''
showinfo=false
zoom=0.9
dotsize=%s
ovalsize=10
groupsize=4
usefoldchange=false
avgfoldchange=false
colorcutoffs=0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;
colorarr=(230;230;230):(207;207;207):(184;184;184):(161;161;161):(138;138;138):(115;115;115):(92;92;92):(69;69;69):(46;46;46):(23;23;23):
</param>""" % (n, args.pvalue, str(dotsize))

    def add_ids(self, ids):
        t = '\n<seq>\n'
        for i in ids:
            t += '>' + i + '\n'
            t += 'X' + '\n'
        t += '</seq>\n'
        if len(ids) != self.n:
            # print 'n != ids'
            raise Exception('n != ids')
        self.txt += t

    def dist_from_matrix(self, rmsds, matrix=0, use_pv=False, use_input_values=False, dont_calc=False, debug=False):
        if dont_calc:
            print('Everything but the dists are generated. Use it to edit the original clans input file.')
            return # for some hardcore debugging ;-)
        t = '\n<hsp>\n'
        c = 0
        c2 = 0
        if debug: print(rmsds)
        for row in rmsds:
            for rmsd in row:
                if c != c2:
                    if use_input_values:
                        dist = float(str(rmsd).replace('e', 'E'))
                    elif use_pv:
                        dist = pv.get_p_value(rmsd, 1 * 38)[0]  # r.get_rmsd_to(r2), 3)
                    else:
                        # 1e-06 10-1 = 9 10-10 0
                        dist = '1.0E-' + str(int(math.floor(matrix.max()) - int(float(rmsd))))
                    t += str(c) + ' ' + str(c2) + ':' + str(dist) + '\n'
                c2 += 1
            c2 = 0
            c += 1

        t += '</hsp>\n'

        if not use_input_values:
            max = math.ceil(matrix.max())
            min = matrix[matrix>0].min()

            self.comment = '# max: %f min (non-zero): %f\n' % (max, min)
            self.comment += '# 1/4 ' + str((max - min) / 4) + ' ' + str(round((max - min) / 4, 0)) + '\n'
            self.comment += '# 1/2 ' + str((max - min) / 2) + ' ' + str(round((max - min) / 2, 0)) + '\n'
            self.comment += '# 1/3 ' + str(((max - min) / 4 ) * 3 ) + ' ' + str(round(((max - min) / 4) * 3, 0)) + '\n'
            for i in range(1,20):
                self.comment += '# connected points with RMSD lower than %iA 1.0E-%i\n' % (i, math.ceil(matrix.max()) - i)
            # 1E-11 = 0
            # 1E-10 = 1-0
            # 1E-9 = 2-1
        self.txt += t
        print(t)

    def dist_from_matrix_mp(self, output_pmatrix_fn, max, min, lines, pmat=False, use_pv=False, use_input_values=False, debug=False):
        if debug:
            print('Everything but the dists are generated. Use it to edit the original clans input file.')
            return # for some hardcore debugging ;-)
        t = '\n<hsp>\n'
        myp = ''
        c = 0
        c2 = 0
        for l in lines:
            for rmsd in l:
                if c != c2:
                    if use_input_values:
                        dist = rmsd.replace('e', 'E')
                    if use_pv:
                        dist = pv.get_p_value(rmsd, 1 * 38)[0]  # r.get_rmsd_to(r2), 3)
                    else:
                        # 1e-06 10-1 = 9 10-10 0
                        dist = '1.0E-' + str(int(math.floor(matrix.max()) - int(float(rmsd))))
                    t += str(c) + ' ' + str(c2) + ':' + str(dist) + '\n'
                    myp += ' ' + str(dist)
                else:
                    myp += ' ' + '0.0'
                c2 += 1
            myp += '\n'
            c2 = 0
            c += 1

        t += '</hsp>\n'

        if not use_input_values:
            max = math.ceil(matrix.max())
            min = matrix[matrix>0].min()

            self.comment = '# max: %f min (non-zero): %f\n' % (max, min)
            self.comment += '# 1/4 ' + str((max - min) / 4) + ' ' + str(round((max - min) / 4, 0)) + '\n'
            self.comment += '# 1/2 ' + str((max - min) / 2) + ' ' + str(round((max - min) / 2, 0)) + '\n'
            self.comment += '# 1/3 ' + str(((max - min) / 4 ) * 3 ) + ' ' + str(round(((max - min) / 4) * 3, 0)) + '\n'
            for i in range(1,20):
                self.comment += '# connected points with RMSD lower than %iA 1.0E-%i\n' % (i, math.ceil(matrix.max()) - i)
            # 1E-11 = 0
            # 1E-10 = 1-0
            # 1E-9 = 2-1

        self.txt += t

        if pmat:
            with open(output_pmatrix_fn, 'w') as f:
                f.write(myp)
        return t

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    """
    https://stackoverflow.com/questions/42908334/checking-if-a-matrix-is-symmetric-in-numpy
    """
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrixfn', help="matrix")

    parser.add_argument('--groups-auto', help="define into how many groups make automatically", type=int)

    parser.add_argument('--color-by-homolog', help="color the same homolog in the same way", action='store_true')
    parser.add_argument('--one-target', help="color the same homolog in the same way", action='store_true')

    parser.add_argument('--shape-by-source', help="shape points based on source, SimRNA vs Rosetta (Rosetta models have 'min' in name')",
                        action='store_true')

    parser.add_argument('--debug', action='store_true')

    parser.add_argument('--groups-dot-size', type=int, default=8)

    parser.add_argument('--dont-calc', action='store_true', help="A simple and dirty trick to get "
                        "generates everything but the main distances from the matrix"
                        "useful if you want to run the script to generate different settings, such as"
                        "colors, groups etc. Run the script and then replace parts of the original "
                        "file with the matrix")

    parser.add_argument('--groups', help="groups, at the moment up to 7 groups can be handle"
                        "--groups 1:native+100:zmp+100:zaa+100:zbaa+100:zcp+100:znc"
                        "--groups 1:native+100:nats+100:hom1+100:hom2+100:hom3+100:hom4"
                        "native will light green"
                        "zmp will be forest green"
                        "(easy to be changed in the future)")

    parser.add_argument('--use-pvalue', action='store_true', help="")
    parser.add_argument('--use-input-values', action='store_true', help="")

    parser.add_argument('--pvalue', default="1.0E-15", help="set p-value for clans.input, default: 1.0E-15")


    parser.add_argument('--output', help="input file for clans, e.g., clans.input", default="clans.input")
    parser.add_argument('--output-pmatrix', action='store_true', help="p value matrix will be saved, see --output-matrix-fn to define name")
    parser.add_argument('--output-pmatrix-fn', default="pmatrix.txt", help="filename of output matrix, pmatrix.txt by default")
    parser.add_argument('--multiprocessing', action='store_true', help="run calculations in parallel way, warning: extra libraries required, see the Docs")

    return parser


#main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    debug = args.debug  # as a short cut later
    f = open(args.matrixfn)
    # OK, check if there is a header = the line with '#'
    headers = f.readline()
    if headers.strip().startswith('#'):
        # if yes, then split remove # and split into lists
        ids = headers.replace('#', '').split()
        for i in ids:
            print(i)
    else:
        # if no, then make a list form [0, # of items in the first line]
        ids = [str(i) for i in range(0, len(headers.split()))]
        f.seek(0)

    # get max
    logging.info(time.strftime("%Y-%m-%d %H:%M:%S"))

    # Collect rmsds into a list
    rmsds = []
    for l in f:
            #for rmsd in map(float,l.split()):  # map(float, s.split())
       		rmsds.append(map(float,l.split()))

    # warning:
    # eh, this matrix is not really used by clanstix main engine
    # it's used to calc min and max value of a matrix
    matrix = np.loadtxt(args.matrixfn)
    if check_symmetric(matrix):
        if args.debug: print('Matrix is symmetrical! ', matrix.shape)
    else:
        raise Exception('Matrix is not symmetrical! Check your matrix', matrix.shape)

    # keep dot quite big by default,
    # but if you use dotsize then keep this one 0
    dotsize =  8
    if args.groups_auto or args.groups:
        dotsize =  0

    c = RNAStructClans(n=len(ids), dotsize=dotsize)  # 200?
    c.add_ids(ids)
    if debug:
        print('dist_from_matrix...')

    if args.multiprocessing:
        import multiprocessing as mp
        #import dill
        #import parmap
        from pathos.multiprocessing import ProcessingPool

        matrix = np.loadtxt(args.matrixfn)
        max = int(math.floor(matrix.max()))
        min = matrix[matrix>0].min()
        clans_list_of_pvalues = ''
        pool = ProcessingPool(mp.cpu_count())
        x=pool.map(c.dist_from_matrix_mp, [args.output_pmatrix_fn], [max], [min], [rmsds], [args.output_pmatrix], )
        pool.close()
        for i in x:
            clans_list_of_pvalues = clans_list_of_pvalues.join(i)
    else:
        c.dist_from_matrix(rmsds, matrix, args.use_pvalue, args.use_input_values, args.dont_calc, args.debug)
        if debug:
            print('process the matrix...')
    #
    # DEFINE GROUPS
    #
    # 1+20+20+20+20+20
    seqgroups = ''
    #colors = ['0;255;102;255', # ligthgreen 1
    #          '0;102;0;255', # forest 2
    # A list of colors used later ....
    colors = [
              '255;102;102;255', # red 3
              '51;51;255;255', # blue 4
              '0;255;255;255', # light blue +1
              '180;38;223;255', # violet 5
              '64;64;64;255', # grey 6
              '255;128;0;255', # orange 7
              '240;230;140;255', #khaki
              '210;105;30;255', # chocolate
              '0;255;255;255',  # cyjan
              '128;0;128;255', # purple 8
              '0;128;128;255', # Teal 9
              '128;0;0;255', # maroon 10

               ]

    # This is pretty much the same list as above, but in here I have more distinguishable colors,
    # as I should be used with less number of groups
    colors_homologs = [
              '255;102;102;255', # red 3
              '51;51;255;255', # blue 4
              '64;64;64;255', # grey 6
              '128;0;128;255', # purple 8
              '128;0;0;255', # maroon 10
              '0;255;255;255',  # cyjan
              '237;41;57;255', # red
              #'210;105;30;255', # chocolate
               ]

    # OK, this is the trick
    # if args.groups_auto is on, then you built args.groups automatically,
    # and then just simply run the next if as it was args.groups in the arguments
    # pretty cool ;-)
    if args.groups_auto:
        # arsg_groups = ''
        from collections import OrderedDict
        # collect groups
        #groups = []
        groups = OrderedDict()
        for i in ids:
            #
            # collect homologs gmp_
            # simrna and farna
            print(i)
            group_name = i[:args.groups_auto]
            if group_name in groups:
                groups[group_name] += 1
            else:
                groups[group_name] = 1
        groups_str = ''
        for g in groups:
            groups_str += str(groups[g]) + ':' + g + '+'
        #
        # this is the trick, you
        args.groups = groups_str[:-1]
        print(groups_str)
        # change this to get 1:hccp+10:zmp+10:xrt

    if args.groups:
        # this is a manual way how to do it
        groups = args.groups.split('+')
        seqgroups = '<seqgroups>\n'
        curr_number = 0
        homologs_colors = {}
        if args.debug: print(args.groups)
        for index, group in enumerate(groups):
            # parse groups
            # type and size will be changed for native
            size = args.groups_dot_size
            dottype = 0
            color = '' # use for diff color selection if tar
            if ':' in group:
                nstruc, name = group.split(':')
                if name == 'native':
                    dottype = 8
                    size = 20
                    color = '0;128;128;255' # olive

                if 'solution' in name:
                    dottype = 8
                    size = 30  # solution is bigger
                    color = '0;128;128;255' # olive

                if args.one_target: # target is target, dont distinguish it
                    if 'tar' in name and 'tar_min' not in name: # SimRNA
                        dottype = 7
                        size = size # 7  # size of SimRNA reference seq points
                        color = '0;255;102;255' # ligthgreen 1

                    if 'tar_min' in name:
                        dottype = 9
                        size = size # 7 # size of Rosetta reference seq points
                        color = '0;102;0;255' # forest 2
                else:
                    if 'tar' in name: # SimRNA
                        size = size #8  # size of SimRNA reference seq points
                        dottype = 0
                        color = '0;255;102;255' # ligthgreen 1

                if args.shape_by_source:
                    # Rosetta models are diamond now
                    if 'min' in name:
                        dottype = 2

                # color by homolog
                if args.color_by_homolog:
                    # 10:gxx_6bd266+10:gxx_min.ou+10:gbaP_d2b57+10:gbaP_min.o+10:gbx_00de79+10:gbx_min.ou+10:gapP_d9d22+10:gapP_min.o+10:gmp_faa97e+10:gmp_min.ou
                    tmp = name.split('_')
                    homolog_name = tmp[0]
                    if debug: 'homolog_name', homolog_name
                    # ignore tar and solution, their colors are defined above ^
                    ## [mm] simrna5x100farna5x100$ git:(master) ✗ rna_clastix.py --groups-auto 8 --color-by-homolog --shape-by-source rp14sub_ref_mapping_refX.txt input2.clans --debug
                    ## 2019-01-09 12:23:21
                    ## 100:tar_min.+100:tar_rp14+1:solution+100:cy2_min.+100:cy2_r14a+100:aj6_min.+100:aj6_r14a
                    ## {'cy2': '210;105;30;255'}
                    ## {'cy2': '210;105;30;255'}
                    ## {'cy2': '210;105;30;255', 'aj6': '0;255;255;255'}
                    ## {'cy2': '210;105;30;255', 'aj6': '0;255;255;255'}
                    ## # max: 18.000000 min (non-zero): 0.810000
                    if homolog_name == 'tar' or homolog_name.startswith('solution'):
                        pass
                    else:
                        if homolog_name in homologs_colors:
                            color = homologs_colors[homolog_name]
                        else:
                            homologs_colors[homolog_name] = colors_homologs.pop()
                            color = homologs_colors[homolog_name]
                        if debug: print(homologs_colors)
            else:
                nstruc = group
                name = 'foo'

            # color is set in the args.color_by_homologs
            # craft seqgroups
            seqgroups += '\n'
            seqgroups += "name=%s\n" % name

            # color hack
            if color:  # this color is fixed for native right now
                seqgroups += "color=%s\n" % color
            else:
                # if beyond index, use different shape
                #try:
                #print(index, len(colors))
                if index >= len(colors):
                    index = index - len(colors) # reset color index
                    dottype = 6  # set dottype when the colors are done 
                seqgroups += "color=%s\n" % colors[index]

            seqgroups += "type=%i\n" % dottype
            # color hack
            seqgroups += "size=%i\n" % size
            seqgroups += "hide=0\n"

            if debug:
                print('name: ' + name + ' ' + colors[index])
            # get numbers - convert nstruc into numbers in Clans format 0;1; etc.
            # remember: it starts from 0
            # --groups 1:hccp+10:zmp+10:xrt
            # hccp
            # [0]
            # zmp
            # [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            # xrt
            # [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            numbers = range(curr_number, curr_number + int(nstruc))  # 0, 1
            curr_number = curr_number + int(nstruc)
            seqgroups += "numbers=%s;\n" % ';'.join([str(number) for number in numbers])
        seqgroups += '</seqgroups>\n'

    with open(args.output, 'w') as f:
        f.write(c.txt)
        if args.multiprocessing:
            f.write(clans_list_of_pvalues)
        f.write(seqgroups)
        f.write(c.comment)
    # if debug: print(c.txt)
    print(c.comment)

    logging.info(time.strftime("%Y-%m-%d %H:%M:%S"))
    if debug:
        print(seqgroups)
