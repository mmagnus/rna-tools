#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Python parser to 3dna <http://x3dna.org/>.

Installation::

  # install the code from http://forum.x3dna.org/downloads/3dna-download/
  Create a copy of the rna_x3dna_config_local_sample.py (remove "_sample") present in rna-tools/rna_tools/tools/rna_x3dna folder.
  Edit this line :
  BINARY_PATH = <path to your x3dna-dssr file>
  matching the path with the path of your x3dna-dssr file.
  e.g. in my case: BINARY_PATH = ~/bin/x3dna-dssr.bin

For one structure you can run this script as::

    [mm] py3dna$ git:(master) ✗ ./rna_x3dna.py test_data/1xjr.pdb
    test_data/1xjr.pdb
    >1xjr nts=47 [1xjr] -- secondary structure derived by DSSR
    gGAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU
    ..(((((((...((((.((((.....))..))..))).).)))))))

For multiple structures in the folder, run the script like this::

    [mm] py3dna$ git:(master) ✗ ./rna_x3dna.py test_data/*
    test_data/1xjr.pdb
    >1xjr nts=47 [1xjr] -- secondary structure derived by DSSR
    gGAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU
    ..(((((((...((((.((((.....))..))..))).).)))))))
    test_data/6TNA.pdb
    >6TNA nts=76 [6TNA] -- secondary structure derived by DSSR
    GCGGAUUUAgCUCAGuuGGGAGAGCgCCAGAcUgAAgAPcUGGAGgUCcUGUGtPCGaUCCACAGAAUUCGCACCA
    (((((((..((((.....[..)))).((((.........)))).....(((((..]....))))))))))))....
    test_data/rp2_bujnicki_1_rpr.pdb
    >rp2_bujnicki_1_rpr nts=100 [rp2_bujnicki_1_rpr] -- secondary structure derived by DSSR
    CCGGAGGAACUACUG&CCGGCAGCCU&CCGGAGGAACUACUG&CCGGCAGCCU&CCGGAGGAACUACUG&CCGGCAGCCU&CCGGAGGAACUACUG&CCGGCAGCCU
    [[[[(((.....(((&{{{{))))))&(((((((.....(.(&]]]]).))))&[[[[[[......[[[&))))]]].]]&}}}}(((.....(((&]]]]))))))

"""
import re
import argparse

from subprocess import Popen, PIPE
from os import remove, path
from rna_tools.rna_tools_config import X3DNA, X3DNA_FP
# x3dna-dssr-64bit


class x3DNAMissingFile(Exception):
    pass


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', '--compact',  action='store_true')
    parser.add_argument('-l', '--show-log',  action='store_true', help="show full log")
    parser.add_argument('files', help='file', nargs='+')
    return parser


class x3DNA(object):

    """
    Atributes:

     **curr_fn**
     **report**

    """

    def __init__(self, pdbfn, show_log=False):
        """Set self.curr_fn based on pdbfn"""
        self.curr_fn = pdbfn
        self.run_x3dna(show_log)
        self.clean_up()

    def __get_report(self):
        """Off right now. Run find_pair

        ..warning:: To get modification use get_modifications()
        """
        cmd = X3DNA_FP + ' ' + self.curr_fn + ' stdout | analyze stdin'  # hack
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = out.stdout.read()
        outerr = out.stderr.read()

        text = ''
        for l in outerr.split('\n'):
            if l.startswith('Time used') or l.startswith('..') or l.startswith('handling') or l.startswith('uncommon residue'):
                continue
            if l.strip():
                text += l + '\n'
        return text.strip()

    def get_modifications(self):
        """Run find_pair to find modifications.
        """
        cmd = X3DNA_FP + ' -p ' + self.curr_fn + ' /tmp/fpout'  # hack
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        outerr = out.stderr.read()

        text = ''
        for l in outerr.split('\n'):
            if l.startswith('uncommon residue'):
                text += l.replace('uncommon residue ', '') + '\n'
        return text.strip()

    def run_x3dna(self, show_log=False):
        """
        """
        cmd = X3DNA + ' -i=' + self.curr_fn
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = str(out.stdout.read().decode())
        outerr = str(out.stderr.read().decode())

        f = open('py3dna.log', 'w')
        f.write(cmd + '\n' + stdout)

        if show_log:
            print(stdout)
        f.close()

        if outerr.find('does not exist!') > -1:  # not very pretty
            raise x3DNAMissingFile
        if outerr.find('not found') > -1:  # not very pretty
            raise Exception('x3dna not found!')

        rx = re.compile('no. of DNA/RNA chains:\s+(?P<no_DNARNAchains>\d+)\s').search(stdout)
        if rx:
            no_of_DNARNA_chains = int(rx.group('no_DNARNAchains'))
            msg = 'py3dna::no of DNARNA chains'
            self.report = msg + '\n' + msg + '\n'  # hack!
        else:
            raise Exception('no_of_DNARNA_chains not found')

        if no_of_DNARNA_chains:
            self.report = stdout.strip()

    def get_ion_water_report(self):
        """@todo
File name: /tmp/tmp0pdNHS
    no. of DNA/RNA chains: 0 []
    no. of nucleotides:    174
    no. of waters:         793
    no. of metals:         33 [Na=29, Mg=1, K=3]
        """
        pass

    def clean_up(self, verbose=False):
        files_to_remove = [
            'dssr-helices.pdb',
            'dssr-pairs.pdb',
            'dssr-torsions.dat',
            'dssr-hairpins.pdb',
            'dssr-multiplets.pdb',
            'dssr-stems.pdb',
            'dssr-Aminors.pdb',
            'hel_regions.pdb',
            'bp_order.dat',
            'bestpairs.pdb',
        ]

        for f in files_to_remove:
            try:
                remove(f)
                pass
            except OSError:
                if verbose:
                    print('can not remove %s' % f)

    def get_seq(self):
        """Get sequence.

        Somehow 1bzt_1 x3dna	UCAGACUUUUAAPCUGA, what is P?
        P -> u
        """
        return self.report.split('\n')[-2].replace('P', 'u').replace('I', 'a')

    def get_secstruc(self):
        """Get secondary structure.
        """
        # return open('dssr-2ndstrs.dbn').read().strip()
        return self.report.split('\n')[-1]


# name
if __name__ == '__main__':
    if not X3DNA:
        raise Exception(
            'Set up BINARY_PATH in rna_x3dna_config_local.py, .e.g "/Users/magnus/work/opt/x3dna/x3dna-dssr"')

    # get parser and arguments
    parser = get_parser()
    args = parser.parse_args()

    # try:
    #    compact = sys.argv[2]
    #    if compact == '-c':
    #       compact = True
    #
    # except IndexError:
    #    compact = False

    for f in args.files:
        if args.compact:
            p = x3DNA(f)
            print((f, p.get_secstruc()))
        else:
            print('# ' + f + ' #')
            p = x3DNA(f, args.show_log)
            # s = p.get_seq()
            # print s
            #s = p.get_secstruc()
            #print(s)
