#!/usr/bin/python
# -*- coding: cp1252 -*-
#
# dihedral_bins.py
#
# Assigns suite names to RNA structures.
#

__author__ = "Kristian Rother, Raphael Bauer"
__credits__ = ["Marcin Domagalski","Magdalena Musielak", "Janusz Bujnicki", "Marie Curie"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

from suite import Suite, SuiteError, SuiteResidueError, \
     SuiteIncompleteResidueError, SuiteBfactorError,\
     SuiteDistanceError
from math import pi,cos
from Scientific.Geometry import *
from dihedral import dihedral, dihedral_from_vectors
from PDB import PDBParser
import os,sys
from suitename import PDBProcessor

BFACTOR_LIMIT = 60
CHECK_BFACTOR = False
CHECK_DISTANCES = False


class EtaThetaError(Exception): pass

class EtaThetaName(Suite):

   def __init__(self):
       Suite.__init__(self)
       self.message = ''

       # suite properties       
       self.eta = None
       self.theta = None
       self.assigned = None # assigned cluster

   def __str__(self):
      result = ''
      if self.assigned:
         result += "%s  assigned\n"%(self.assigned.name)
      result += Suite.__str__(self) + '\n'
      result += self.message +'\n'
      return result
                     
      # assign to one of 12 delta,delta,gamma bins
      self.ddg_bin = "%s%s%s"%(self.puckerdm,self.puckerd,self.gammaname)
      self.message += "ddg==%s. "%self.ddg_bin

   def set_residues(self, resi1, resi2, resi3):
       """
       Takes two PDB.Residue objects and calculates suite dihedrals
       from them. If the residues do not apply, a SuiteError is raised.
       """
       rv1 = self.get_vectors(resi1, ["C4*","P"])
       rv2 = self.get_vectors(resi2, ["C4*","P"])
       rv3 = self.get_vectors(resi3, ["P"])
       self.resi1 = resi1
       self.resi2 = resi2
       self.resi3 = resi3
       self.vecs1 = rv1
       self.vecs2 = rv2
       self.vecs3 = rv3
       # calculate pseudotorsion angles
       self.eta  = dihedral_from_vectors(rv1["C4*"],rv2["P"],rv2["C4*"],rv3["P"])
       self.theta = dihedral_from_vectors(rv1["P"],rv1["C4*"],rv2["P"],rv2["C4*"])

   def get_et_bin(self):
      """
      Assigns characters depending on eta/theta dihedral angles, and returns
      a string of two characters, where
      A -  0.. 9.999°
      B - 10..19.999°
      C - 20... etc.
      """
      ebin = 65+int(self.eta)//10
      tbin = 65+int(self.theta)//10
      result = chr(ebin) + chr(tbin)
      return result
   
      
   def assign_name(self,check_bfactor=CHECK_BFACTOR,check_distances=CHECK_DISTANCES):
      """
      Returns the name of the suite, if it can be assigned.
      """
      try:
         # self.confirm(False,check_bfactor, check_distances)
         self.assigned = self.get_et_bin()
      except EtaThetaError,e:
         self.message += str(e)
         raise e


class EtaThetaNameProcessor(PDBProcessor):

   def report(self):
         keys = self.counts.keys()
         keys.sort()
         for k in keys:
            print k,self.counts[k]
         print "suites assigned",snp.good
         print "bad residues   ",snp.resierrors

   def get_string_for_chain(self,chain):
      result = ''
      first = None
      second = None
      third = None
      for resi in chain.child_list:
         first, second, third = second, third, resi
         if first and second and third:
            try:
               et = EtaThetaName()
               et.set_residues(first, second, third)
               # phosphate atom of third residue
               rv3 = et.get_vectors(third, ['P'])
               et.assign_name(rv3['P'])
               self.good += 1
               result += et.assigned
            except SuiteIncompleteResidueError, e:
               self.resierrors += 1
               result += '--'            
      return result

if __name__ == '__main__':
   print """
erena - Suite Suffix Tree
(c) 2008 Kristian Rother and Raphael Bauer
create strings from eta-theta angle binning.

usage: dihedral_bins.py <path with pdb files> <outfile>
   """   
   snp = EtaThetaNameProcessor()
   path = sys.argv[1]
   outfile = sys.argv[2]
   snp.loop_path(path)
   snp.report()
   snp.write_output(outfile)
   
   

