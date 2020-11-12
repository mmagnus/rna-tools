#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Bio.PDB.PDBParser
import Bio.PDB.Superimposer
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO, Superimposer

import argparse

import optparse
import sys
import math
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
        #self.atoms = []
        #if save:
        #    self.save() # @save

    def __get_atoms(self):
        self.atoms=[]
        # [<Atom P>, <Atom C5'>, <Atom O5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C2'>, <Atom O2'>, <Atom C1'>, <Atom N1>, <Atom C2>, <Atom N3>, <Atom C4>, <Atom C5>, <Atom C6>, <Atom N6>, <Atom N7>, <Atom C8>, <Atom N9>, <Atom OP1>, <Atom OP2>]
        for res in self.struc.get_residues():
            for at in res:
                #print(res.resname, res.get_id, at)
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

    def get_rmsd_to(self, other_rnamodel, way="", save=True):
        """Calc rmsd P-atom based rmsd to other rna model"""
        sup = Bio.PDB.Superimposer()

        if way == 'backbone+sugar':
            self.atoms_for_rmsd = []
            for a in self.atoms:
                if a.name in "P,OP1,OP2,C5',O5',C4',O4',C3',O3',C2',O2',C1'".split(','):
                    self.atoms_for_rmsd.append(a)

            other_atoms_for_rmsd = []
            for a in other_rnamodel.atoms:
                if a.name in "P,OP1,OP2,C5',O5',C4',O4',C3',O3',C2',O2',C1'".split(','):
                    other_atoms_for_rmsd.append(a)
        else:
            self.atoms_for_rmsd = self.atoms
            other_atoms_for_rmsd = other_rnamodel.atoms

        sup.set_atoms(self.atoms_for_rmsd, other_atoms_for_rmsd)
        rms = round(sup.rms, 3)

        if save:
            ## io = Bio.PDB.PDBIO()
            ## sup.apply(self.struc.get_atoms())
            ## io.set_structure(self.struc)
            ## io.save("aligned.pdb")

            io = Bio.PDB.PDBIO()
            sup.apply(other_rnamodel.struc.get_atoms())
            io.set_structure(other_rnamodel.struc)
            io.save(other_rnamodel.fpath.replace('.pdb', '_' + args.suffix + '.pdb'))
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
    parser.add_argument('--way', help="e.g., backbone+sugar")
    parser.add_argument("-s", "--save", action="store_true", help="set suffix with --suffix, by default: aligned")
    parser.add_argument('files', help='files', nargs='+')
    parser.add_argument("--debug", action="store_true")
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    target = RNAmodel(args.target)

    models = args.files  # opts.input_dir
    tmp = []
    if args.ignore_files:
        for f in args.files:
            if args.debug: print(f)
            if args.ignore_files in f:
                continue
            tmp.append(f)
        models = tmp
    print('# of models:', len(models))
    c = 1
    t = 'model,rmsd\n'
    for m in models:
        mrna = RNAmodel(m)
        #print r1.fn, r2.fn, r1.get_rmsd_to(r2)#, 'tmp.pdb')
        rmsd = target.get_rmsd_to(mrna, args.way, args.save) #, 'tmp.pdb')
        #print rmsd
        t += mrna.fn + ',' + str(rmsd) + '\n'
        #break    
    print(t.strip())
    if args.result:
        with open(args.result, 'w') as f:
            f.write(t)
