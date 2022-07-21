#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Examples::

  triplexibility.py -t 2nd_triplex_FB_ABC.pdb  triples-all-v2-rpr/*.pdb --way backbone+sugar --save --suffix ABC --result abc.csv

  triplexibility.py -t 2nd_triplex_FB_1AUA3_rpr.pdb test_data/triples2/Triple_cWW_tSH_GCA_exemplar_rpr_ren.pdb  --way backbone+sugar --save --column-name 'rmsd' --triple-mode > out.txt

"""
from __future__ import print_function

import Bio.PDB.PDBParser
import Bio.PDB.Superimposer
import copy
import warnings
warnings.filterwarnings('ignore', '.*Invalid or missing.*',)
warnings.filterwarnings('ignore', '.*with given element *',)
warnings.filterwarnings('ignore', '.*is discontinuous*',)
warnings.filterwarnings('ignore', '.*Ignoring unrecognized record*',)

import sys
from icecream import ic
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')

from rna_tools.rna_tools_lib import set_chain_for_struc, RNAStructure
from rna_tools.rna_tools_config import RT

import argparse
import sys
import glob
import re
import os

class RNAmodel:
    def __init__(self, fpath):
        # parser 1-5 -> 1 2 3 4 5
        self.struc = Bio.PDB.PDBParser().get_structure('', fpath)
        self.__get_atoms()
        self.fpath = fpath
        self.fn = os.path.basename(fpath)
        # self.atoms = []
        # if save:
        #    self.save() # @save

    def __get_atoms(self):
        self.atoms=[]
        # [<Atom P>, <Atom C5'>, <Atom O5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom O2'>, <Atom C1'>, <Atom N1>, <Atom C2>, <Atom N3>, <Atom C4>, <Atom C5>, <Atom C6>, <Atom N6>, <Atom N7>, <Atom C8>, <Atom N9>, <Atom OP1>, <Atom OP2>]
        for res in self.struc.get_residues():
            for at in res:
                #if args.debug: print(res.resname, res.get_id, at)
                self.atoms.append(at)
        return self.atoms

    def __str__(self):
        return self.fn #+ ' # beads' + str(len(self.residues))

    def __repr__(self):
        return self.fn #+ ' # beads' + str(len(self.residues))

    def get_report(self):
        """Str a short report about rna model""" 
        t = ' '.join(['File: ', self.fn, ' # of atoms:', str(len(self.atoms)), '\n'])
        for r,a in zip(self.residues, self.atoms ):
            t += ' '.join(['resi: ', str(r) ,' atom: ', str(a) , '\n' ])
        return t

    def get_rmsd_to(self, other_rnamodel, way="", triple_mode=False, save=True, tseq=''):
        """Calc rmsd P-atom based rmsd to other rna model

        sugar now 10 atoms ;-) """
        sup = Bio.PDB.Superimposer()
        if way in ['c1', 'backbone+sugar', 'sugar']:
            if way == 'c1':
                atomslist = ["C1'"]# ,"C2'","O4'"] #, "C2'"]

            elif way == 'sugar':
                atomslist = "C5',O5',C4',O4',C3',O3',C2',O2',C1'".split(',')

            elif way == 'backbone+sugar':
                atomslist = "P,OP1,OP2,C5',O5',C4',O4',C3',O3',C2',O2',C1'".split(',')

            self.atoms_for_rmsd = []
            for a in self.atoms:
                nt = a.get_parent().get_resname().strip()
                if nt in ['G', 'A']:
                    atomslistx = atomslist + ['N9']
                if nt in ['C', 'U']:
                    atomslistx = atomslist + ['N1']
                if a.name in atomslistx:
                    self.atoms_for_rmsd.append(a)
            if args.debug: print('atoms_for_rmsd', len(self.atoms_for_rmsd))

            other_atoms_for_rmsd = []
            for a in other_rnamodel.atoms:
                nt = a.get_parent().get_resname().strip()
                if nt in ['G', 'A']:
                    atomslistx = atomslist + ['N9']
                if nt in ['C', 'U']:
                    atomslistx = atomslist + ['N1']
                if a.name in atomslistx:
                    other_atoms_for_rmsd.append(a)

            if args.debug: print('other_atoms_for_rmsd', len(other_atoms_for_rmsd))
                                 
        elif way == 'c1+Nx':
            self.atoms_for_rmsd = []
            for a in self.atoms:
                nt = a.get_parent().get_resname().strip() # G
                if nt in ['G', 'A']:
                    atomslist = ["C1'", 'N9'] # , 'N1']
                if nt in ['C', 'U']:
                    atomslist = ["C1'", 'N1'] # , 'N1']
                if a.name in atomslist:
                    self.atoms_for_rmsd.append(a)

            if args.debug: print('atoms_for_rmsd', len(self.atoms_for_rmsd))

            other_atoms_for_rmsd = []
            for a in other_rnamodel.atoms:
                nt = a.get_parent().get_resname().strip() # G
                if nt in ['G', 'A']:
                    atomslist = ["C1'", 'N9'] # , 'N1']
                if nt in ['C', 'U']:
                    atomslist = ["C1'", 'N1'] # , 'N1']
                if a.name in atomslist:
                    other_atoms_for_rmsd.append(a)
            if args.debug: print('other_atoms_for_rmsd', len(other_atoms_for_rmsd))

        else:
            self.atoms_for_rmsd = self.atoms
            other_atoms_for_rmsd = other_rnamodel.atoms

        # calc rmsd #
        if not tseq:
            tseq = ''
            rt = None
            for a in self.atoms_for_rmsd:
                r = a.get_parent()
                if r != rt:
                    rt = r
                    tseq += r.get_resname().strip()

        if triple_mode:
            def chunks(lst, n):
                """Yield successive n-sized chunks from lst.
                https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
                """
                for i in range(0, len(lst), n):
                    yield lst[i:i + n]

            rmsd_min = 10000 # ugly
            import itertools
            # ok, for different residues now it's a problem  
            per = list(itertools.permutations([0, 1, 2])) # [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
            lst = list(chunks(other_atoms_for_rmsd,
                              int(len(other_atoms_for_rmsd)/3))) # for 1 atom, this will be 1 x 3 [3 residues]
                       # so so len is 3 atoms so / by 3 to get how many atoms per residue
            print(lst)

            sup_min = None
            seq_min = 'not yet obtained, rmsd rejected!'
            p_min = None

            rms = -1

            for p in per:
                patoms = []
                for i in p: # p=(1, 2, 3)
                    patoms.extend(lst[i])

                #print(self.atoms_for_rmsd)
                ## print('patoms')
                ## for a in patoms:
                ##     print(a, a.get_parent().get_id())
                ## print('self.atoms_for_rmsd')
                ## for a in self.atoms_for_rmsd:
                ##     print(a, a.get_parent().get_id())
                #sup.set_atoms(patoms, self.atoms_for_rmsd)

                rt = None
                seq = ''
                for a in patoms:
                    r = a.get_parent()
                    if r != rt:
                        rt = r
                        seq += r.get_resname().strip()
             
                ic(tseq.lower(), seq.lower(), other_rnamodel.fpath)
                # dont' even calc rmsd if the curr seq and tseq are not the same
                if tseq.lower() == seq.lower():  # only if seq is the same
                    print(self.atoms_for_rmsd)
                    print(patoms)
                    print(len(self.atoms_for_rmsd), len(patoms))
                    for a, b in zip(self.atoms_for_rmsd, patoms):
                        print(a, a.get_parent().id, a.get_parent().get_resname(),
                           b.get_parent().id, b.get_parent().get_resname())
                    print(len(self.atoms_for_rmsd), len(patoms))
                    try:
                        sup.set_atoms(self.atoms_for_rmsd, patoms)
                    except:
                        pass
                    rms = round(sup.rms, 2)
                    if rms < rmsd_min:
                            rmsd_min = rms
                            sup_min = copy.copy(sup)
                            suffix = seq
                            p_min = p
                            seq_min = seq
                            if args.debug: 'set new rmsd_min', rmsd_min

                    if args.debug:
                        print(p, '', [i + 1 for i in p], end=' ')
                        print(seq, 'seq_min: ' + seq_min, rms)

            if p_min: # found target sequence
                # what is this? ;-)
                index = [0 ,0 ,0]
                index[0] = p_min.index(0)
                index[1] = p_min.index(1)
                index[2] = p_min.index(2)

                # ugly re-set 123 to crazy id ! + 100, so this will
                # fill up 1 2 3 for the second for
                rs = other_rnamodel.struc[0].get_residues()
                for i, r in enumerate(rs):
                    r.id = (' ', index[i] + 253, ' ') # ugly, some random offset

                for i, r in enumerate(other_rnamodel.struc[0].get_residues()):
                    r.id = (' ', index[i] + 1, ' ')
                    if args.debug: print('r', r)

                io = Bio.PDB.PDBIO()
                sup_min.apply(other_rnamodel.struc.get_atoms())
                # if args.debug: print(p_min, [i + 1 for i in p_min])
                io.set_structure(other_rnamodel.struc)

                args.save_here = True
                if args.save_here and save:
                    folder = os.path.basename(self.fpath.replace('.pdb', '_' + args.folder_prefix + '_aligned'))
                    # print(f)
                    try:
                        os.mkdir(folder)
                    except:
                        pass
                    fout = folder + os.sep + "{:1.2f}".format(rmsd_min) + '-' + os.path.basename(other_rnamodel.fpath)#.replace('.pdb', '-' + str(rms) + '.pdb'))
                    #_s' + suffix + '.pdb'))
                else:
                    fout = other_rnamodel.fpath.replace('.pdb', '_aligned.pdb')#_s' + suffix + '.pdb')

                if args.debug: print(fout)

                if save:
                    io.save(fout)
                    # ugly set chain to A
                    set_chain_for_struc(fout, 'A', save_file_inplace=True)
                    # and now run this to sort into 1 2 3

                    r = RNAStructure(fout)
                    remarks = r.get_remarks_text()
                    r1 = r.get_res_text('A', 1)
                    r2 = r.get_res_text('A', 2)
                    r3 = r.get_res_text('A', 3)
                    with open(fout, 'w') as f:
                        f.write(remarks)
                        f.write(r1)
                        f.write(r2)
                        f.write(r3)
                    r.reload()
                    r.get_rnapuzzle_ready()
                    if rmsd_min < 1:  # !!!!!!!!!!!! ugly
                         r.write()
                return str(rmsd_min)# + ',s' + seq_min + ',' + os.path.basename(fout)
        else:
            # check if number of the same atoms #
            # if not the same then return -1    #
            # print(len(self.atoms_for_rmsd), len(other_atoms_for_rmsd))
            if len(self.atoms_for_rmsd) != len(other_atoms_for_rmsd):
                return -1
            sup.set_atoms(self.atoms_for_rmsd, other_atoms_for_rmsd)
            rms = round(sup.rms, 2)

            if save:
                ## io = Bio.PDB.PDBIO()
                ## sup.apply(self.struc.get_atoms())
                ## io.set_structure(self.struc)
                ## io.save("aligned.pdb")

                io = Bio.PDB.PDBIO()
                sup.apply(other_rnamodel.struc.get_atoms())
                io.set_structure(other_rnamodel.struc)

            args.save_here = True
            if args.save_here and save:
                f = os.path.basename(self.fpath.replace('.pdb', '_aligned'))
                # print(f)
                try:
                    os.mkdir(f)
                except:
                    pass
                # fout = f + os.sep + os.path.basename(other_rnamodel.fpath.replace('.pdb', '_aligned_s' + suffix + '.pdb'))
                fout = f + os.sep + os.path.basename(other_rnamodel.fpath)#.replace('.pdb', '_aligned'))
            else:
                fout = other_rnamodel.fpath.replace('.pdb', '_aligned.pdb')
            if save:
               io.save(fout)
        return rms

    
def get_rna_models_from_dir(directory):
    models = []
    if not os.path.exists(directory):
        raise Exception('Dir does not exist! ', directory)
    files = glob.glob(directory + "/*.pdb")
    files_sorted = sort_nicely(files)
    for f in files_sorted:
        #print f
        models.append(RNAmodel(f))
    return models

def sort_nicely( l ):
   """ Sort the given list in the way that humans expect.

   http://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
   """
   convert = lambda text: int(text) if text.isdigit() else text
   alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
   l.sort( key=alphanum_key )
   return l


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t',"--target",
                         help="target file")
    parser.add_argument('--result', #default='out.rmsd',
                         help="result file")
    parser.add_argument('--ignore-files', default='aligned', help="use to ignore files, .e.g. with 'aligned'")
    parser.add_argument('--suffix', default='aligned', help="used with --saved, by default: aligned")
    parser.add_argument('--way', help="e.g., backbone+sugar, c1, c1+Nx")
    parser.add_argument('--triple-mode', help="same crazy strategy to align triples", default=True, action="store_true")
    parser.add_argument('--column-name', help="name column for rmsd, by default 'rmsd', but can also be a name of the target file", default="rmsd")
    parser.add_argument("-s", "--save", action="store_true", help="set suffix with --suffix, by default: aligned")
    parser.add_argument("--folder-prefix", default = '', help="folder name, t2-3-UAU_NAME_aligned")

    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--sort", action="store_true", help='sort results based on rmsd (ascending), --result must be set up')
    parser.add_argument("--tseq", help='target sequence, e.g. acu, find only triples of this sequence [use the order for the seq taken from the input PDB file, literally order of residues in a pdb file]')
    parser.add_argument('--files', help='files', nargs='+', default=RT + '/rna_tools/tools/triplexibility/db/triples-all-v2-rpr/*.pdb')
    return parser

# main
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    if not args.debug:
        ic.disable()

    target = RNAmodel(args.target)

    # a trick to get files be default if there is only a path (so a string ;-))
    if not isinstance(args.files, str):
        models = args.files
    else:
        import glob
        models = glob.glob(args.files) # opts.input_dir

    ## tmp = []
    ## if args.ignore_files:
    ##     for f in args.files:
    ##         if args.debug: print(f)
    ##         if args.ignore_files in f:
    ##             continue
    ##         tmp.append(f)
    ##     models = tmp
 
    print('# of models:', len(models))

    c = 1
    t = 'fn,' + args.column_name + '\n' # ,aligned_seq, aligned_fn\n'
    for m in models:
        mrna = RNAmodel(m)
        #print r1.fn, r2.fn, r1.get_rmsd_to(r2)#, 'tmp.pdb')
        # print(m)
        rmsd = target.get_rmsd_to(mrna, args.way, args.triple_mode, args.save, args.tseq) #, 'tmp.pdb')
        #except:
        if 0:
            print(m)
            sys.exit(1)
        #print rmsd
        t += mrna.fn + ',' + str(rmsd) + '\n'
        #break    

    print(t.strip())

    if args.result:
        with open(args.result, 'w') as f:
            f.write(t)

        if args.sort:
            import pandas as pd
            df = pd.read_csv(args.result)
            df = df[df.rmsd != -1]
            df = df.sort_values('rmsd')
            df.to_csv(args.result, index=False)
            print('saved: %s' % args.result)
            print(df.to_string(index=False))
