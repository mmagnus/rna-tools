#!/usr/bin/python
# -*- coding: cp1252 -*-
#
# suite.py
#
# Stores suites from RNA structures.
#

__author__ = "Kristian Rother, Raphael Bauer"
__credits__ = ["Marcin Domagalski","Magdalena Musielak", "Janusz Bujnicki", "Marie Curie"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"


from rna_tools.tools.mini_moderna3.moderna.PDB.PDBParser import PDBParser
from rna_tools.tools.mini_moderna3.moderna.dihedral import dihedral

FIRST_SUITEATOMS = ["C5*","C4*","C3*","O3*"]
SECOND_SUITEATOMS = ['P',"O5*","C5*","C4*","C3*","O3*"]

class SuiteError(Exception): pass
class SuiteResidueError(SuiteError): pass
class SuiteIncompleteResidueError(SuiteResidueError): pass
class SuiteBfactorError(SuiteResidueError): pass
class SuiteDistanceError(SuiteResidueError): pass


class Suite(SuiteAngles):
   """
   Class that reads dihedral angles out of Bio.PDB.Residue objects.
   Does the validation for suites.
   """
   def __init__(self):
       SuiteAngles.__init__(self,EMPTY_ANGLES)
       self.resi1 = None
       self.resi2 = None
       self.vecs1 = None
       self.vecs2 = None

       self.puckerdm = '' # pucker of 5' half-suite
       self.puckerd = ''  # pucker of 3' half-suite 
       self.gammaname = ''

   def set_residues(self, resi1, resi2):
       """
       Takes two PDB.Residue objects and calculates suite dihedrals
       from them. If the residues do not apply, a SuiteError is raised.
       """
       rv1 = self.get_vectors(resi1, FIRST_SUITEATOMS)
       rv2 = self.get_vectors(resi2, SECOND_SUITEATOMS)
       self.resi1, self.resi2 = resi1, resi2
       self.vecs1, self.vecs2 = rv1, rv2
       # calculate torsion angles
       deltam  = dihedral(rv1["C5*"],rv1["C4*"],rv1["C3*"],rv1["O3*"])
       epsilon = dihedral(rv1["C4*"],rv1["C3*"],rv1["O3*"],rv2["P"])
       zeta    = dihedral(rv1["C3*"],rv1["O3*"],rv2["P"],rv2["O5*"])
       alpha   = dihedral(rv1["O3*"],rv2["P"],rv2["O5*"],rv2["C5*"])
       beta    = dihedral(rv2["P"],rv2["O5*"],rv2["C5*"],rv2["C4*"])
       gamma   = dihedral(rv2["O5*"],rv2["C5*"],rv2["C4*"],rv2["C3*"])
       delta   = dihedral(rv2["C5*"],rv2["C4*"],rv2["C3*"],rv2["O3*"])
       angles = [deltam, epsilon, zeta, alpha, beta, gamma, delta]
       SuiteAngles.set_angles(self,angles)

   def confirm_distances(self):
      """
      Raises an exception if any interatomic distances are above 2.2.
      """
      # build list of vectors in right order
      vectors = [self.vecs1[name] for name in FIRST_SUITEATOMS]
      vectors += [self.vecs2[name] for name in SECOND_SUITEATOMS]
      # check distances
      for i in range(len(vectors)-1):
         dist = vectors[i+1] - vectors[i]
         if dist.norm() > 2.2:
            raise SuiteDistanceError("interatomic distance of %5.2f too big"%dist.norm())
            
   def confirm_bfactors(self):
      """
      Raises a SuiteBfactorError if any of the atoms in
      the list of atom names has a bfactor too high.
      """
      rn = ((self.resi1,FIRST_SUITEATOMS),(self.resi2,SECOND_SUITEATOMS))
      for resi, atom_names in rn:
         for name in atom_names:
            if resi[name].bfactor >= BFACTOR_LIMIT:
               raise SuiteBfactorError('too big B-factor: %6.3f'%resi[name].bfactor)

   def confirm(self,check_angles,check_bfactors,check_distances):
       if check_angles: self.confirm_angles()
       if check_bfactors: self.confirm_bfactors()
       if check_distances: self.confirm_distances()

   def get_pucker(self, delta):
      if ALLOWED_ANGLES['delta3min'] <= delta <= ALLOWED_ANGLES['delta3max']:
         return '3'
      elif ALLOWED_ANGLES['delta2min'] <= delta <= ALLOWED_ANGLES['delta2max']:
         return '2'

   def filter_epsilon(self):
      """Checks whether epsilon is in the allowed range."""
      if self.epsilon < ALLOWED_ANGLES['epsilonmin'] or self.epsilon > ALLOWED_ANGLES['epsilonmax']:
         raise SuiteTriageError("epsilon outlier. ")

   def filter_deltam(self):
      """Checks whether delta-1 is in the allowed range."""
      self.puckerdm = self.get_pucker(self.deltam)
      if not self.puckerdm:
        raise SuiteTriageError("bad deltam. ")

   def filter_delta(self):
      """Checks whether delta is in the allowed range."""
      self.puckerd = self.get_pucker(self.delta)
      if not self.puckerd:
        raise SuiteTriageError("bad delta. ")

   def filter_gamma(self):
      """Checks whether gamma is in the allowed range."""
      if(ALLOWED_ANGLES['gammapmin'] <= self.gamma <= ALLOWED_ANGLES['gammapmax']):
           self.gammaname = 'p'
      elif(ALLOWED_ANGLES['gammatmin'] <= self.gamma <= ALLOWED_ANGLES['gammatmax']):
           self.gammaname = 't'
      elif(ALLOWED_ANGLES['gammammin'] <= self.gamma <= ALLOWED_ANGLES['gammammax']):
           self.gammaname = 'm'
      else:
           raise SuiteTriageError("gamma outlier. ")

   def filter_alpha(self):
      """Checks whether alpha is in the allowed range."""
      if(ALLOWED_ANGLES['alphamin'] > self.alpha or self.alpha > ALLOWED_ANGLES['alphamax']):
         raise SuiteTriageError("alpha outlier. ")

   def filter_beta(self):
      """Checks whether beta is in the allowed range."""
      if(ALLOWED_ANGLES['betamin'] > self.beta or self.beta > ALLOWED_ANGLES['betamax']):
         raise SuiteTriageError("beta outlier. ")

   def filter_zeta(self):
      """Checks whether zeta is in the allowed range."""
      if(ALLOWED_ANGLES['zetamin'] > self.zeta or self.zeta > ALLOWED_ANGLES['zetamax']):
         raise SuiteTriageError("zeta outlier. ")

   def filter_angles(self):
      self.filter_epsilon()
      self.filter_deltam()
      self.filter_delta()
      self.filter_gamma()
      self.filter_alpha()
      self.filter_beta()
      self.filter_zeta()
    
   def get_vectors(self,resi, atom_names):
      """Creates a dict of Vector objects from a residue object."""
      vectors = {}
      for name in atom_names:
         if resi.has_id(name):
            vectors[name] = resi[name].get_vector()
         else:
            raise SuiteIncompleteResidueError('atom %s does not exist.'%name)
      return vectors


