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

.. warning:: This script should not be used in this given form with Parallel because it process output files from x3dna that are named always in the same way, e.g. dssr-torsions.txt. #TODO

"""
import re
import argparse

from subprocess import Popen, PIPE
import os
from rna_tools.rna_tools_config import X3DNA, X3DNA_FP
from icecream import ic;import sys;ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True,contextAbsPath=True);ic.configureOutput(prefix='')
# x3dna-dssr-64bit


class x3DNAMissingFile(Exception):
    pass


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', '--compact',  action='store_true')
    parser.add_argument('-x', '--rerun',  action='store_true')
    parser.add_argument('-s', '--show',  action='store_true', help="show results")
    parser.add_argument('--pymol',  action='store_true', help='get resi to color code puckers in PyMOL')
    parser.add_argument('-l', '--show-log',  action='store_true', help="show full log")
    parser.add_argument('-v', '--verbose',  action='store_true', help="show full log")
    parser.add_argument('-f', '--foldability',  type=int, default=50, help="move files with foldability < this to folder _low_foldability_X")
    parser.add_argument('files', help='file', nargs='+')
    return parser


class x3DNA(object):

    """
    Atributes:

     **curr_fn**
     **report**

    """

    def __init__(self, pdbfn, show_log=False, verbose=False):
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

    def run_x3dna(self, show_log=False, verbose=False):
        """
        """
        cmd = X3DNA + ' -i=' + self.curr_fn
        out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)

        stdout = str(out.stdout.read().decode())
        outerr = str(out.stderr.read().decode())

        # print(stdout, outerr, cmd)
        #f = open('py3dna.log', 'w')
        #if verbose: print(f'cmd: {cmd}')
        #f.write(cmd + '\n' + stdout)
        # if show_log:
        #     print(stdout)
        # f.close()

        if outerr.find('does not exist!') > -1:  # not very pretty
            raise x3DNAMissingFile
        if outerr.find('not found') > -1:  # not very pretty
            raise Exception(f'x3dna not found! {self.curr_fn} {cmd}')

        if verbose:
            print(stdout)
        self.report = stdout # msg + '\n' + msg + '\n'  # hack!
        return        

        rx = re.compile('no. of DNA/RNA chains:\s+(?P<no_DNARNAchains>\d+)\s').search(stdout)
        if rx or True:
            no_of_DNARNA_chains = int(rx.group('no_DNARNAchains'))
            msg = 'py3dna::no of DNARNA chains'
            self.report = msg + '\n' + msg + '\n'  # hack!
        else:
            print(stdout)
            raise Exception(f'no_of_DNARNA_chains not found  {self.curr_fn} {cmd}')

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
                os.remove(f)
                pass
            except OSError:
                if verbose:
                    print('error: can not remove %s' % f)

    def get_seq(self):
        """Get sequence.

        Somehow 1bzt_1 x3dna	UCAGACUUUUAAPCUGA, what is P?
        P -> u
        
        """
        #return self.report.split('\n')[-2].replace('P', 'u').replace('I', 'a')
        pass

    def get_secstruc(self):
        """Get secondary structure.
        """
        def percent_paired(dot_bracket):
            # Remove invalid characters (e.g. '&')
            valid_chars = [c for c in dot_bracket if c in '.()']

            total = len(valid_chars)
            ic(total, valid_chars.count('(') + valid_chars.count(')'))
            paired = valid_chars.count('(') + valid_chars.count(')')  # each base counts, not each pair
            if total == 0:
                return 0.0  # avoid division by zero

            return round(100 * paired / total, 2)

        def parse_chain_block(lines):
            for i in range(len(lines) - 2):
                if "[chain]" in lines[i]:
                    header = lines[i].strip()
                    seq = lines[i + 1].strip()
                    struct = lines[i + 2].strip()
                    return header, seq, struct
            return None, None, None

        header, seq, struct = parse_chain_block(self.report.split('\n'))
        #print("Header:", header)
        #print("Sequence:", seq)
        #print("Structure:", struct)
        if struct:
            return struct, percent_paired(struct)
        else:
            return struct, None
        hits = re.search("as a whole and per chain.*?\n(?P<ss>.+?)\n\*", self.report, re.DOTALL|re.MULTILINE)
        if hits:
             return hits.group('ss').strip()
        else:
            ic(self.report)
            return self.report.split('\n')[-1] # tofix

    def get_torsions(self, outfn) -> str:
        """Get torsion angles into 'torsion.csv' file::

            nt,id,res,alpha,beta,gamma,delta,epsilon,zeta,e-z,chi,phase-angle,sugar-type,ssZp,Dp,splay,bpseq
            1,g,A.GTP1,nan,nan,142.1,89.5,-131.0,-78.3,-53(BI),-178.2(anti),358.6(C2'-exo),~C3'-endo,4.68,4.68,29.98,0
            2,G,A.G2,-75.8,-167.0,57.2,79.5,-143.4,-69.7,-74(BI),-169.2(anti),5.8(C3'-endo),~C3'-endo,4.68,4.76,25.61,0
        
        """
        angles = ''
        save = False

        c2pendo = []
        c3pendo = []
        if not os.path.isfile('dssr-torsions.txt'):
            
            print(f'Problem to get torsion angles for {self.curr_fn}')
            return f'Problem to get torsion angles for {self.curr_fn}'
            
        for l in open('dssr-torsions.txt'):
            if 'nt               alpha    beta' in l:
                save = True
                l = l.replace('nt',  'nt id   res   ')
            if '***************' in l and save:
                save = False
            if save:
                if "~C2'-endo" in l:
                    c2pendo.append(l.split()[0])
                if "~C3'-endo" in l:
                    c3pendo.append(l.split()[0])
                angles += l

        c2 = 'color pink, resi ' + '+'.join(c2pendo)
        c3 = 'color blue, resi ' + '+'.join(c3pendo)

        if args.pymol:
            print(c2)
            print(c3)
            return
            
        import re
        nangles = ''
        #'9 C    41',
        #'10 C     0',
        bpseq = ['bpseq'] + [x.strip().split()[2] for x in open('dssr-2ndstrs.bpseq').readlines()]
        for i, l in enumerate(angles.split('\n')):
            if l.strip():
                l = re.sub(r'\s+', ',', l, 0, re.MULTILINE)
                if bpseq[i] == '0':
                    bpseq[i] = 'no paired'
                else:
                    bpseq[i] = 'paired'
                l = l[1:] + ',' + bpseq[i] + '\n'
                nangles += l
        nangles = re.sub(r'---', 'nan', nangles, 0, re.MULTILINE)
        with open(outfn, 'w') as f:
            f.write(nangles.strip())

        if args.show:
            import pandas as pd
            df = pd.read_csv(outfn)
            print(df)
        return nangles.strip()
    
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
            #print(f'input: {f}')
            outfn = os.path.basename(f.replace('.pdb', '')) + '-torsion-paired.csv'
            if not args.rerun:
                if os.path.isfile(outfn):
                    print(f'skip: {f}, use --rerun to run analysis again')
                    continue
            p = x3DNA(f, args.show_log, args.verbose)
            s = p.get_seq()
            #print(s)
            s, foldability = p.get_secstruc()
            print(f'{f} {s} {foldability}%')
            #s = p.get_torsions(outfn)
            #if args.verbose: print(s)
            import math
            def round_up_to_10(x: int) -> int:
                return int(math.ceil(x / 10.0)) * 10

            p.clean_up(args.verbose)
            if foldability is not None and foldability < args.foldability:
                print(f'Warning: {f} has low foldability {foldability}%')
                # foldibilty very 5And move to folder
                foldability = round_up_to_10(foldability)
                os.makedirs(f'_low_foldability_{foldability}', exist_ok=True)
                os.rename(f, os.path.join(f'_low_foldability_{foldability}', os.path.basename(f)))
                
                
                
