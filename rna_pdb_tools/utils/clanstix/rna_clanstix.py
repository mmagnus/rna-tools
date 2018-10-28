#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Clanstix - a tool for visualizing RNA 3D structures based on pairwise structural similarity with Clans.

We hacked Clans thus instead of BLAST-based distances between sequences, you can analyze distances between structures described as p-values of rmsd (based on the method from the Dokholyan lab.)

Running Clans:
To run CLANS you need to have Java 1.4 or better installed (java can be downloaded HERE). For full functionality you will also need the NCBI BLAST,PSI-BLAST and formatdb executables (NCBI). For command line parameters and basic help please refer to the README file.
(source: http://www.eb.tuebingen.mpg.de/research/departments/protein-evolution/software/clans.html)

.. image:: ../../rna_pdb_tools/utils/clanstix/doc/yndSrLTb7l.gif

The RMSDs between structures are converted into p-values based on the method from the Dokholyan lab.

Color groups
---------------------------------------
You can color your groups:

.. image:: ../../rna_pdb_tools/utils/clanstix/doc/rna_clanstix.png

To get colors, run a cmd like this::

   rna_clastix.py rnapz17_matrix_farfar_HelSeedCst.txt --groups 20:seq1+20+20+20+20+20+20:seq10

where with the ``+`` sign you separate groups. Each group has to have a number of structures. Optionally it can have a name, e.g., ``20:seq1``, use ``:`` as a separator. If a provided name is ``native`` then this group will be shown as starts.

Get inspiration for more colors (http://www.rapidtables.com/web/color/RGB_Color.htm)

How to use ClanstixRNA?
----------------------------------------

1. Get a matrix of distances, save it as e.g. matrix.txt
2. run ClanstixRNA on this matrix to get an input file to Clans (e.g. clans_rna.txt)::

     clanstix.py test_data/matrix.txt > clans_run.txt

3. open CLANS and click File -> Load run and load clans_run.txt
4. You're done! :-)

Hajdin, C. E., Ding, F., Dokholyan, N. V, & Weeks, K. M. (2010). On the significance of an RNA tertiary structure prediction. RNA (New York, N.Y.), 16(7), 1340–9. doi:10.1261/rna.1837410

An output of this tool can be viewed using CLANS.

Frickey, T., & Lupas, A. (2004). CLANS: a Java application for visualizing protein families based on pairwise similarity. Bioinformatics (Oxford, England), 20(18), 3702–4. doi:10.1093/bioinformatics/bth444
"""

import argparse
import rna_pdb_tools.utils.rmsd_signif.rnastruc_pred_signif as pv
import numpy as np
import math


class RNAStructClans:
    """

    Usage::

        >>> f = open('matrix.txt')
        >>> ids = f.readline().replace('#','').split()
        >>> c = RNAStructClans(n=len(ids)) # 200?
        >>> c.add_ids(ids)
        >>> c.dist_from_matrix(f)
        >>> print(c.txt)
    """

    def __init__(self, n=10):
        self.n = n
        self.comment = ''
        self.txt = """sequences=%i
<param>
maxmove=0.1
pval=1.0E-15
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
cluster2d=false
blastpath=blastall -p blastp -I T -T T
formatdbpath=/home/sdh/apps/blast-2.2.17/bin/formatdb
showinfo=false
zoom=1.0
dotsize=1
ovalsize=10
groupsize=4
usefoldchange=false
avgfoldchange=false
colorcutoffs=0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;
colorarr=(230;230;230):(207;207;207):(184;184;184):(161;161;161):(138;138;138):(115;115;115):(92;92;92):(69;69;69):(46;46;46):(23;23;23):
</param>""" % n

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

    def dist_from_matrix(self, lines, matrix=0, use_pv=False):
        t = '\n<hsp>\n'
        c = 0
        c2 = 0
        for l in lines:
            for rmsd in l.split():
                if c != c2:
                    if use_pv:
                        dist = pv.get_p_value(rmsd, 1 * 38)[0]  # r.get_rmsd_to(r2), 3)
                    else:
                        # 1e-06 10-1 = 9 10-10 0
                        dist = '1.0E-' + str(int(math.floor(matrix.max()) - int(float(rmsd))))
                    t += str(c) + ' ' + str(c2) + ':' + str(dist) + '\n'
                c2 += 1
            c2 = 0
            c += 1

        t += '</hsp>\n'

        self.comment = '# max: %f min (non-zero): %f\n' % (math.ceil(matrix.max()), matrix[matrix>0].min())
        for i in range(1,20):
            self.comment += '# connected points with RMSD lower than %iA 1.0E-%i\n' % (i, math.ceil(matrix.max()) - i)
        t += '# max: %f min (non-zero): %f' % (math.ceil(matrix.max()), matrix[matrix>0].min())
        t += ' # 1A RMSD range is for lower than 1.0E-%f' % (math.ceil(matrix.max()) - 1)
        t += ' # 2A RMSD range is for lower than 1.0E-%f' % (math.ceil(matrix.max()) - 2)
        # 1E-11 = 0
        # 1E-10 = 1-0
        # 1E-9 = 2-1
        self.txt += t


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrixfn', help="matrix")
    parser.add_argument('--groups', help="groups, at the moment up to 7 groups can be handle (easy to be changed in the future)")
    parser.add_argument('--pvalue', action='store_true',
                        help="")
    return parser


# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    f = open(args.matrixfn)
    ids = f.readline().replace('#', '').split()
    # print ids
    # print len(ids)
    # get max
    matrix = np.loadtxt(args.matrixfn)

    c = RNAStructClans(n=len(ids))  # 200?
    c.add_ids(ids)
    c.dist_from_matrix(f, matrix, args.pvalue)
    print(c.txt)

    # print color
    # 1+20+20+20+20+20
    seqgroups = ''
    colors = ['0;255;102;255', # ligthgreen
              '255;102;102;255', # red
              '0;102;0;255', # forest
              '51;51;255;255', # blue
              '180;38;223;255', # violet
              '64;64;64;255', # grey
              '255;128;0;255' # orange
                  ]
    if args.groups:
        groups = args.groups.split('+')
        seqgroups = '<seqgroups>\n'
        curr_number = 0
        for index, group in enumerate(groups):
            # parse groups
            # type and size will be changed for native
            size = 10
            dottype = 0
            if ':' in group:
                nstruc, name = group.split(':')
                if name == 'native':
                    dottype = 8
                    size = 20
            else:
                nstruc = group
                name = 'foo'

            # craft seqgroups
            seqgroups += "name=%s\n" % name
            seqgroups += "type=%i\n" % dottype
            seqgroups += "color=%s\n" % colors[index]
            seqgroups += "size=%i\n" % size
            seqgroups += "hide=0\n"
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
    print(seqgroups)
    print(c.comment)
