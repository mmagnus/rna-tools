#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""rna_refinement - RNA refinement with QRNAS.

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

Be default the script searches QRNAS in <rna-tools>/opt/qrnas/ .

Usage of QRNA::

    QRNA - Quick Refinement of Nucleic Acids (0.2 alpha)
         by Juliusz Stasiewicz (jstasiewicz@genesilico.pl)

    To use type:
      QRNA -i <input PDBfile> [-o <output PDBfile>] [-c <configfile>] [-p] [-m <restraintsfile>]
    OR specify <input PDBfile>, <output PDBfile> and <restraintsfile> in <configfile> and type just:
      QRNA -c <configfile>

**Installation of this util**

Set up in your bashrc::

   export QRNAS_PATH=<your path to qrnas> # e.g. /home/magnus/src/rna-tools/opt/qrnas

but default rna-tools searches for qrnas in <rna-tools>/opt/qrnas.

**QRNAS at Peyote2**

There is no problem to run QRNAS at our Genesilico cluster, `peyote2`. Tested by mmagnus --170822. Copy files of QRNAS to peyote and run ``./qrnamake sequential``.

To run it at a cluster with the Sun Grid Engine queuing system (this one with qusb ;-))::

     for p in *.pdb; do echo "rna_refinement.py $p >& ${p}.log" | qsub -cwd -V -pe mpi 1 -N "r_$p" ; done

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
from shutil import copyfile
from rna_tools.rna_tools_config import QRNAS_PATH, QRNAS_CONFIG_PATH

class QRNAS:
    """QRNAS"""
    def run(self, inputfile, outputfile, run_in_tmp=False, job_id_random=False, steps=10, interactive=False, verbose=False):
        """Run QRNAS.

        Args:
           inputfile     (str) : path to a input file, use .pdb extensions
           outputfile    (str) : path to on output file
           run_in_tmp    (bool): if yes, run in /tmp otherwise run in currect-directory/tmp/THEjHxilN3nLx2Aj8REg/input
           job_id_random (bool): if yes, then job id will be like, e.g. tmp/gOIFfSdo9tnelFtvs3A7/output.pdb, if now
                                 tmp/output/output.pdb; output, not input because then you can run the same input
                                 a view times and get different outputs (and name of outputs will give
                                 folder names
           steps         (int) : # of steps
           interactive   (bool): if yes, use os.system (see progress for cmd on the screen) or subprocess (wait till it's finished)

        Returns:
           none: works on input/output files

        """
        cwd = os.getcwd()
        # get config
        print(QRNAS_PATH)

        #QRNAS_FF_DIR
        
        conftxt = open(QRNAS_PATH + os.sep + 'configfile.txt').read()
        conftxt_tmp = re.sub('\#?\s?NSTEPS.+\d+', 'NSTEPS   ' + str(steps), conftxt)
        conftxt_tmp = re.sub('\#?\s?VERBOSE.+\d+', 'VERBOSE 1', conftxt)

        if job_id_random:
            JOB_ID = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(20))
        else:
            JOB_ID = os.path.basename(outputfile.replace('.pdb', ''))

        if run_in_tmp:
            JOB_PATH = '/tmp/' + os.sep + JOB_ID + os.sep
        else:
            JOB_PATH = cwd + os.sep + 'tmp' + os.sep + JOB_ID + os.sep  # run it in place?

        try:
            os.makedirs(JOB_PATH)
        except:
            pass

        # get temp config
        with open(JOB_PATH + os.sep + 'configfile.txt','w') as f:
            f.write(conftxt_tmp)

        # copy input to qrnas folder
        qrnas_inputfile = JOB_PATH + os.path.basename(inputfile)
        qrnas_outputfile = JOB_PATH + os.path.basename(inputfile).replace('.pdb', '_refx.pdb')
        copyfile(inputfile, qrnas_inputfile)

        os.chdir(QRNAS_PATH)
        cmd = 'export QRNAS_FF_DIR="' + QRNAS_PATH + '/forcefield" && ' + QRNAS_PATH + '/QRNA -i ' + qrnas_inputfile + \
              ' -c ' + JOB_PATH + 'configfile.txt ' + \
              ' -o ' + qrnas_outputfile
        if verbose:
            print(cmd)

        stdout = open(JOB_PATH + 'stdout.txt', 'w')
        stderr = open(JOB_PATH + 'stderr.txt', 'w')

        if interactive:
            os.system(cmd)
        else:
            subprocess.call(cmd, shell=True, stdout=stdout, stderr=stderr)
            if verbose:
                with open(JOB_PATH + 'stderr.txt') as f:
                    print(f.read())

        os.chdir(cwd)
        print ("Save to %s" % outputfile)
        copyfile(qrnas_outputfile, outputfile)


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--steps', help="# of steps, default: 20k ", default=20000)
    parser.add_argument('fn', help="input pdb file")
    parser.add_argument('-o', '--output_file', help="output pdb file")
    parser.add_argument('-i', '--interactive', help="", default=False,
                              action="store_true")
    parser.add_argument('-v', '--verbose', help="",
                              action="store_true")
    return parser

# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    q = QRNAS()
    if not args.output_file:
        output_file = args.fn.replace('.pdb', '_refx.pdb')
    else:
        output_file = args.output_file
    q.run(args.fn, output_file, steps=args.steps, interactive=args.interactive, verbose=args.verbose)
