#!/usr/bin/python

import Bio.PDB.PDBParser
import Bio.PDB.Superimposer
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO, Superimposer

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
        for res in self.struc.get_residues():
                self.atoms.extend(res.get_list())
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

    def get_rmsd_to(self, other_rnamodel, output=''):
        """Calc rmsd P-atom based rmsd to other rna model"""
        sup = Bio.PDB.Superimposer()
        sup.set_atoms(self.atoms, other_rnamodel.atoms)
        rms = round(sup.rms, 3)

        if output:
            io = Bio.PDB.PDBIO()
            sup.apply(self.struc.get_atoms())
            io.set_structure( self.struc )
            io.save("aligned.pdb")

            io = Bio.PDB.PDBIO()
            sup.apply(other_rnamodel.struc.get_atoms())
            io.set_structure( other_rnamodel.struc )
            io.save("aligned2.pdb")
            
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

if __name__ == '__main__':
    optparser=optparse.OptionParser(usage="%prog [<options>]")

    optparser.add_option('-i',"--input_dir", type="string",
                         dest="input_dir",
                         default='',
                         help="")

    optparser.add_option('-o',"--matrix_fn", type="string",
                         dest="matrix_fn",
                         default='matrix.txt',
                         help="ouput, matrix")

    optparser.add_option("-s", "--save",
                     action="store_true", default=False, dest="save", help="")

    
    (opts, args)=optparser.parse_args()

    if len(sys.argv) == 1:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)

    input_dir = opts.input_dir
    matrix_fn = opts.matrix_fn

    models = get_rna_models_from_dir(input_dir)        

    print '# of models:', len(models)

    f = open(matrix_fn, 'w')
    t = '# '
    for r1 in models:
        #print r1,
        t += str(r1) + ' '
    #print
    t += '\n'

    c = 1
    for r1 in models:
            for r2 in models:
                #print
                #print r1.fn, r2.fn, r1.get_rmsd_to(r2)#, 'tmp.pdb')
                rmsd = r1.get_rmsd_to(r2) #, 'tmp.pdb')
                #print rmsd
                t += str(rmsd) + ' '
                #break    
            #print
            t += '\n'
            
    f.write(t)
    f.close()

    print t.strip() # matrix

    if True:
        print 'matrix was created! ', matrix_fn
    else:
        print 'matrix NOT was created!'
