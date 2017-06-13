#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Clanstix - a tool for visualizing RNA 3D structures based on pairwise structural similarity with Clans.

We hacked Clans thus instead of BLAST-based distances between sequences, you can analyze distances between structures described as p-values of rmsd (based on the method from the Dokholyan lab.)

Running Clans:
To run CLANS you need to have Java 1.4 or better installed (java can be downloaded HERE). For full functionality you will also need the NCBI BLAST,PSI-BLAST and formatdb executables (NCBI). For command line parameters and basic help please refer to the README file.
(source: http://www.eb.tuebingen.mpg.de/research/departments/protein-evolution/software/clans.html)

.. image:: ../../rna_pdb_tools/utils/clanstix/doc/yndSrLTb7l.gif

The RMSDs between structures are converted into p-values based on the method from the Dokholyan lab.

How to use ClanstixRNA?

1. Get a matrix of distances, save it as e.g. matrix.txt
2. run ClanstixRNA on this matrix to get an input file to Clans (e.g. clans_rna.txt)::

     clanstix.py test_data/matrix.txt > clans_run.txt

3. open CLANS and click File -> Load run and load clans_run.txt
4. You're done! :-)

Hajdin, C. E., Ding, F., Dokholyan, N. V, & Weeks, K. M. (2010). On the significance of an RNA tertiary structure prediction. RNA (New York, N.Y.), 16(7), 1340–9. doi:10.1261/rna.1837410

An output of this tool can be viewed using CLANS.

Frickey, T., & Lupas, A. (2004). CLANS: a Java application for visualizing protein families based on pairwise similarity. Bioinformatics (Oxford, England), 20(18), 3702–4. doi:10.1093/bioinformatics/bth444
"""

import random
import sys
import rna_pdb_tools.utils.rmsd_signif.rnastruc_pred_signif as pv

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
        self.txt = """sequences=%i
<param>
maxmove=0.1
pval=1.0E-8
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
showinfo=true
zoom=1.0
dotsize=10
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
            #print 'n != ids'
            raise Exception('n != ids')
        self.txt += t

    def dist_from_matrix(self, lines):
        t = '\n<hsp>\n'
        c = 0
        c2 = 0
        for l in lines:
            for i in l.split():
                if c != c2:
                    #t += str(c) + ' ' + str(c2) + ':' + i.upper() + '\n' #'1E-' + str(random.randint(0,10)) + '\n' # #i.upper() + '\n'
                    #if random.randint(0,2):
                    rms = i
                    dist = pv.get_p_value(rms, 1 * 38)[0] #r.get_rmsd_to(r2), 3)
                    #print i, dist
                    t += str(c) + ' ' + str(c2) + ':' + str(dist) + '\n' # '1E-' + str(random.randint(0,15)) + '\n' # #i.upper() + '\n'
                c2 += 1
            c2 = 0
            c += 1

        t += '</hsp>\n'
        self.txt += t

if __name__ == '__main__':
    f = open(sys.argv[1]) # matrix.txt
    ids = f.readline().replace('#','').split()
    #print ids
    #print len(ids)
    c = RNAStructClans(n=len(ids)) # 200?
    c.add_ids(ids)
    c.dist_from_matrix(f)
    print(c.txt)
