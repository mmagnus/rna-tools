#!/usr/bin/python
# -*- coding: utf-8 -*-

"""RNA refinement with QRNAS

Models of RNA 3D structures obtained by modeling methods often suffer from local inaccuracies such as clashes or physically improbable bond lengths, backbone conformations, or sugar puckers. To ensure high quality of models, a procedure of re nement should be applied as a  nal step in the modeling pipeline. The software tool QRNAS was developed in our laboratory to perform local re nement of nucleic acid structures based on an extended version of the AMBER force field. The extensions consist of energy terms associated with introduction of explicit hydrogen bonds, idealization of base pair planarity and regularization of backbone conformation.

Read more: Piatkowski, P., Kasprzak, J. M., Kumar, D., Magnus, M., Chojnowski, G., & Bujnicki, J. M. (2016). RNA 3D Structure Modeling by Combination of Template-Based Method ModeRNA, Template-Free Folding with SimRNA, and Refinement with QRNAS. Methods in Molecular Biology (Clifton, N.J.), 1490(Suppl), 217-235. http://doi.org/10.1007/978-1-4939-6433-8_14

Right now, there is 20k steps of refinement.

.. image:: ../pngs/qrnas_0k.png

The initial structure, ``179c48aa-c0d3-4bd6-8e06-12081da22998_ALL_thrs6.20A_clust01-000001_AA.pdb``.

.. image:: ../pngs/qrnas_3k.png

after 3k, ~10min

.. image:: ../pngs/qrnas_10k.png

after 10k steps, around 30min

.. image:: ../pngs/qrnas_20k.png

after 20k steps, around 1h.

**Installation of QRNAS**

Download the QRNAS package from http://genesilico.pl/qrnas/,
unzip the archive, and compile it with the following command::

    ./qrnamake sequential

This should create an executable version of QRNAS.

.. warning:: Please, change the name of the binary file from QRNA to QRNAS!

Be default the script searches QRNAS in <rna-pdb-tools>/opt/qrnas/ .

Usage of QRNA::

    QRNA - Quick Refinement of Nucleic Acids (0.2 alpha)
         by Juliusz Stasiewicz (jstasiewicz@genesilico.pl)

    To use type:
      QRNA -i <input PDBfile> [-o <output PDBfile>] [-c <configfile>] [-p] [-m <restraintsfile>]
    OR specify <input PDBfile>, <output PDBfile> and <restraintsfile> in <configfile> and type just:
      QRNA -c <configfile>

**Installation of this util**

Set up in your bashrc::

   export QRNAS_PATH=<your path to qrnas> # e.g. /home/magnus/src/rna-pdb-tools/opt/qrnas

but default rna-pdb-tools searches for qrnas in <rna-pdb-tools>/opt/qrnas.

**QRNAS at Peyote2**

There is no problem to run QRNAS at our Genesilico cluster, `peyote2`. Tested by mmagnus --170822. Copy files of QRNAS to peyote and run ``./qrnamake sequential``.

DONE:

- [x] clean up the output structure
- [x] configuration should not be hardcoded
"""

from __future__ import print_function
import argparse
import re
import os
import subprocess
import random
import string
import sys

from shutil import copyfile

try:
    PATH = os.environ['RNA_PDB_TOOLS']
except:
    print ('Set up RNA_PDB_TOOLS, see Installation note')
    pass
else:
    QRNAS_PATH = os.getenv('QRNAS_PATH', PATH + '/opt/qrnas/')

class QRNAS:
    """QRNAS"""
    def run(self, inputfile, outputfile, steps = 10):
        """
        :param inputfile: 
        :param outputfile: 
        :param steps: 
        """
        cwd = os.getcwd()
        # get config
        conftxt = open(QRNAS_PATH + os.sep + 'configfile.txt').read()
        conftxt_tmp = re.sub('\#?\s?NSTEPS.+\d+', 'NSTEPS   ' + str(steps), conftxt)
        JOB_ID = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(20))

        JOB_PATH = '/tmp/' + os.sep + JOB_ID + os.sep
        os.makedirs(JOB_PATH)

        # get temp config
        with open(JOB_PATH + os.sep + 'configfile.txt','w') as f:
            f.write(conftxt_tmp)

        # copy input to qrnas folder
        qrnas_inputfile = JOB_PATH + os.path.basename(inputfile)
        qrnas_outputfile = JOB_PATH + os.path.basename(inputfile).replace('.pdb', '.refx.pdb')
        copyfile(inputfile, qrnas_inputfile)

        os.chdir(QRNAS_PATH)
        cmd = './QRNAS -i ' + qrnas_inputfile + \
          ' -c ' + JOB_PATH + 'configfile.txt ' + \
          ' -o ' + qrnas_outputfile
        print(cmd)
        subprocess.call(cmd, shell=True)

        os.chdir(cwd)
        print ("Save to %s" % outputfile)
        copyfile(qrnas_outputfile, outputfile)

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--steps', help="# of steps, default: 20k ", default=20000)         
    parser.add_argument('fn', help="input pdb file")     
    parser.add_argument('-o', '--output_file', help="output pdb file")
    return parser
#main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    q = QRNAS()
    if not args.output_file:
        output_file = args.fn.replace('.pdb', '_refx.pdb')
    else:
        output_file = args.output_file
    q.run(args.fn, output_file, steps = args.steps)
