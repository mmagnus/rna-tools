#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""rna_rosetta_run.py - prepare & run ROSETTA simulations

Based on C. Y. Cheng, F. C. Chou, and R. Das, Modeling complex RNA tertiary folds with Rosetta, 1st ed., vol. 553. Elsevier Inc., 2015.
http: // www.sciencedirect.com / science / article / pii / S0076687914000524

The script makes(1) a folder for you job, with seq.fa, ss.fa, input file is copied as input.fa to the folder(2) make helices(3) prepare rosetta input files(4) sends jobs to the cluster.

The header is take from the fast file(`` > /header / ``) not from the filename of your Fasta file.

I discovered this::

     qstat -xml | tr '\n' ' ' | sed 's#<job_list[^>]*>#\n#g' \
>   | sed 's#<[^>]*>##g' | grep " " | column -t

(https://stackoverflow.com/questions/26104116/qstat-and-long-job-names) so there is now need to shorted my job ids.

Helix
-------------------------------------------------------

Run::

    rna_rosetta_run.py -i -e -r -g -c 200 cp20.fa

`-i`::

    # prepare a folder for a run
    >cp20
    AUUAUCAAGAAUCUCAAAGAGAGAUAGCAACCUGCAAUAACGAGCAAGGUGCUAAAAUAGAUAAGCCAAAUUCAAUUGGAAAAAAUGUUAA
    .(((((....(((((.....)))))(((..(((((..[[[[..)).))).)))......))))).((((......)))).......]]]].

    [peyote2] ~ rna_rosetta_run.py -i cp20.fa
    run rosetta for:
    cp20
    AUUAUCAAGAAUCUCAAAGAGAGAUAGCAACCUGCAAUAACGAGCAAGGUGCUAAAAUAGAUAAGCCAAAUUCAAUUGGAAAAAAUGUUAA
    .(((((....(((((.....)))))(((..(((((..[[[[..)).))).)))......))))).((((......)))).......]]]].
    /home / magnus // cp20 / created
    Seq & ss created

Troubleshooting.

If one of the helices is missing you will get::

    IOError: [Errno 2] No such file or directory: 'helix1.out'
    rosetta_submit.py README_FARFAR o 500 100 taf
    Could not find:  README_FARFAR

and the problem was a1 and g8 pairing::

    outputting command line to:  helix0.RUN  # previous helix #0
    Sequence:  AUGG CCGG
    Secstruc:  (((())))
    Not complementary at positions a1 and g8!
    Sequence:  GUGGG CCCAU
    Secstruc:  ((((()))))

    Writing to fasta file:  helix2.fasta  # next helix #2

My case with a modeling of rp12

    Sequence:  cc gc
    Secstruc:  (())
    Not complementary at positions 1 and 4!

edit the secondary structure, run the program with -i(init, to overwrite seq.fa, ss.fa) and then it works.


Notes::

   rp17hc 6 characters

"""
from __future__ import print_function
import argparse
import textwrap
import os
import glob
import subprocess
import shutil
import math

import os
import sys

try:
    from rna_tools.rna_tools_config import RNA_ROSETTA_RUN_ROOT_DIR_MODELING
except:
    RNA_ROSETTA_RUN_ROOT_DIR_MODELING = ''
    print ('Set up rna_rosetta_run_root_dir_for_modeling in rpt_config_local.py')


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=CustomFormatter)

    parser.add_argument('-i', '--init', help='prepare _folder with seq and ss',
                        action='store_true')
    parser.add_argument('-e', '--helices', help='produce h(E)lices',
                        action='store_true')
    parser.add_argument('-r', '--rosetta', help='prepare rosetta files (still you need `go` to send jobs to a cluster)',
                        action='store_true')
    parser.add_argument('-g', '--go', help=go.__doc__,
                        action='store_true')
    parser.add_argument('-m', '--motif', help="include a motif file, e.g. -s E-loop_1q9a_mutated_no_flanks_renumber.pdb")
    parser.add_argument('-n', '--nstruc', help="# of structures you want to get",
                        default=10000, type=int)

    parser.add_argument('-c', '--cpus', help='# of cpus to be used', default=200,
                        type=int)

    parser.add_argument('--sandbox', help="where to run it (default: RNA_ROSETTA_RUN_ROOT_DIR_MODELING",
                        default=RNA_ROSETTA_RUN_ROOT_DIR_MODELING)

    parser.add_argument('file', help=textwrap.dedent(
        """file: \n > a04\nUAUAACAUAUAAUUUUGACAAUAUGGGUCAUAAGUUUCUACCGGAAUACCGUAAAUAUUCUGACUAUGUAUA\n((((.((((...((.((((.......)))).))........(.(((((.......))))).)..))))))))"""))
    return parser


def prepare_folder(args, header, seq, ss, path):
    """Make folder for you job, with seq.fa, ss.fa, input file is copied as input.fa to the folder.

    For ss lowercase is needed when used with motifs, otherwise::

        [peyote2] aa20$ rna_rosetta_run.py -r -m E-loop_1q9a_mutated_no_flanks_renumber_for_acy20.pdb ~/aa20.fa
        2019-05-03 21:31:30,842 rpt_config_local.py::<module>::rpt_config_local loading...
        run rosetta for:
        aa20
        UACGUUCAUCAUCCGUUUGGAUGACGGAAGUAAGCGAAAGCUGAAGGAACGCAUG
        ..(((((.((((((....))))))..[.....(((....)))....)))))]...
        rna_denovo_setup.py -fasta seq.fa -secstruct_file ss.fa -cycles 20000 -no_minimize -nstruct 50  -s E-loop_1q9a_mutated_no_flanks_renumber_for_acy20.pdb  -silent helix0.out helix1.out helix2.out -input_silent_res 3-7 47-51 9-14 19-24 33-35 40-42

        Sequence:  UACGUUCAUCAUCCGUUUGGAUGACGGAAGUAAGCGAAAGCUGAAGGAACGCAUG
        Secstruc:  ..(((((.((((((....))))))..[.....(((....)))....)))))]...
        aaguagaag
        AAGUAGAAG
        Traceback (most recent call last):
          File "/home/magnus/opt/rosetta_src_2016.13.58602_bundle/tools/rna_tools/bin//rna_denovo_setup.py", line 170, in <module>
            raise ValueError('The sequence in %s does not match input sequence!!' % pdb)
        ValueError: The sequence in E-loop_1q9a_mutated_no_flanks_renumber_for_acy20.pdb does not match input sequence!!
        rosetta_submit.py README_FARFAR o 200 100 aa20_
        Could not find:  README_FARFAR
    """
    d = path
    try:
        os.mkdir(d)
        print(d, 'created')
    except OSError:
        print(d, 'created is already created')
        pass

    with open(d + "seq.fa", "w") as f:
        f.write(header + '\n')
        f.write(seq)

    with open(d + "ss.fa", "w") as f:
        f.write(ss.lower())
    print('Seq & ss created')
    shutil.copyfile(args.file, d + 'input.fa')


def prepare_helices():
    """Make helices(wrapper around 'helix_preassemble_setup.py')

    .. warning:: I think multiprocessing of helixX.run does not work."""
    # run helix_p..
    cmd = 'helix_preassemble_setup.py -secstruct ss.fa -fasta seq.fa'
    os.system(cmd)

    # find all helix
    helix_runs = glob.glob('*RUN')
    print(helix_runs)

    f = open('HRUNS', 'w')
    for h in helix_runs:
        f.write(open(h).read().strip() + ' & \n')
    f.close()

    # does not work (!)
    # os.system('chmod +x CMDLINES')
    # os.system('./CMDLINES')
    # ./CMDLINES: 2: source: not found

    # os.system runs multiprocessing, but does not wait for the rest of the program
    # hmm... it does not wait because I used & ???
    # this works in multirpocessing mode but it does not wait for `-r` !!! so if your run -e only it's OK.
    # don't combile -e with -r because making helices will not wait to run -r (!) and  you will get an error
    # and only helices made and then the program will finish
    if False:
        os.system('bash HRUNS')
    else:
        p = subprocess.Popen(open('HRUNS').read(),
                             shell=True, stderr=subprocess.PIPE)
        p.wait()
        stderr = p.stderr.read().strip()
        if stderr:
            print(stderr)


def prepare_rosetta(header, cpus, motif, nstruc):
    """Prepare ROSETTA using rna_denovo_setup.py

    cpus is used to calc nstruc per job to get 10k structures per full run::

    Args:

      nstruc(int): how many structures you want to obtain
      nstruct = int(math.floor(20000 / cpus))
      motif: motif file; e.g., -s E-loop_1q9a_mutated_no_flanks_renumber.pdb
      50 (nstruc) = 10k / 200 (cpus)

    """
    # get list line
    helices = open('CMDLINES').readlines()[-1].replace('#', '')
    njobs = cpus  # 500
    nstruct = int(math.floor(nstruc / cpus))  # 20000/500 -> 40
    if motif:
        cmd_motif = ' -s ' + motif
    else:
        cmd_motif = ''
    cmd = 'rna_denovo_setup.py -fasta seq.fa -secstruct_file ss.fa -cycles 20000 -no_minimize -nstruct ' + \
        str(nstruct) + ' ' + cmd_motif + ' ' + helices
    print(cmd)
    os.system(cmd)
    # change to 50 per job (!)
    # 50 * 100 = 10k ?
    # dont change this 100 (!) it might not run on peyote2 with values like 99999 !
    # cmd is
    # rna_tools/bin//rosetta_submit.py <text file with rosetta command> <outdir> <# jobs> <# hours> <job name>
    # manually: [peyote2] a22$ rosetta_submit.py README_FARFAR o 200 100 a22_
    cmd = 'rosetta_submit.py README_FARFAR o ' + \
        str(njobs) + ' 100 ' + header + '_'
    print(cmd)
    os.system(cmd)


def go():
    """send jobs to a cluster(run qsubMINI)"""
    os.system('chmod +x ./qsubMINI')
    os.system('./qsubMINI')


def main():
    """Pipeline for modeling RNA"""
    args = get_parser().parse_args()

    if args.file:
        f = open(args.file)
        header = f.readline().replace('>', '').replace(' ', '').strip()
        seq = f.readline().strip()
        ss = f.readline().strip()
        cpus = int(args.cpus)
        print('run rosetta for:')
        print(header)

        print(seq)
        print(ss)

        if RNA_ROSETTA_RUN_ROOT_DIR_MODELING.strip() == '':
            print('Set RNA_ROSETTA_RUN_ROOT_DIR_MODELING in your rpt_config file.')
            return

        path = args.sandbox + os.sep + header + \
            os.sep  # RNA_ROSETTA_RUN_ROOT_DIR_MODELING
        curr = os.getcwd()
        if args.init:
            prepare_folder(args, header, seq, ss, path)
        try:
            os.chdir(path)
        except OSError:
            print('You have to make a folder first! use --init')
            sys.exit(1)

        if args.helices:
            prepare_helices()
        if args.rosetta:
            prepare_rosetta(header, cpus, args.motif, args.nstruc)
        if args.go:
            go()

        os.chdir(curr)


# main
if __name__ == '__main__':
    main()
