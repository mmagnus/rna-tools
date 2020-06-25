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


from PDB.PDBParser import PDBParser
from math import pi,cos
from dihedral import dihedral_from_vectors
import os


FIRST_SUITEATOMS = ["C5*","C4*","C3*","O3*"]
SECOND_SUITEATOMS = ['P',"O5*","C5*","C4*","C3*","O3*"]

EMPTY_ANGLES = [-1.0]*7

POWER = 3.0

class SuiteError(Exception): pass
class SuiteResidueError(SuiteError): pass
class SuiteIncompleteResidueError(SuiteResidueError): pass
class SuiteBfactorError(SuiteResidueError): pass
class SuiteDistanceError(SuiteResidueError): pass

class SuiteAngles:
    """
    A class that stores the seven consecutive dihedral angles
    a suite consists of and allows basic calculations with them.
    """
    def __init__(self, angles):
        """
        angles is a list of seven floats.
        [deltam epsilon zeta   alpha  beta   gamma  delta]
        """
        self.set_angles(angles)

    def __getitem__(self, i):
        return self.angles[i]

    def check_positive(self,angle):
        if angle<0: angle += 360.0
        return angle
    
    def set_angles(self, angles):
        """angles is a list of seven floats."""
        angles = map(self.check_positive,angles)
        self.deltam, self.epsilon,self.zeta,self.alpha,\
            self.beta,self.gamma,self.delta = angles
        self.angles = angles

    def get_7d_diff(self, suite_angles):
        """
        Returns SuiteAngleslist with the differences
        between both suites angles.
        """
        angles = [suite_angles[i]-self[i] for i in range(7)]
        return SuiteAngles(angles)

    def dotproduct(self, suite_angles, nang):
        """
        Returns dot product of both suite angles.
        If nang=4 is set, only epsilon, zeta, alpha, beta will be used.
        Otherwise all seven
        """
        indices = nang==4 and range(1,5) or range(7)
        products = [self[i]*suite_angles[i] for i in indices]
        return sum(products)

    def get_hyperellipsoid_dist(self, suite_angles, nang, weights):
       """
       Calculates distance to hyperellipsoid.
       suite contains the center point of the ellipsoid,
       weights the ellisoid radii.
       nange can be 4 or 7 to choose distance in 4D or 7D.
       Returns float
       -> Piet Hein superellipse, hyperellipsoids. 
           0      1      2      3      4      5      6  
        deltam epsilon zeta   alpha  beta   gamma  delta 
           X                                  X      X   
           X not used in 4 angle distance calc
       """
       indices = (nang==4) and range(1,5) or range(7)
       # normalize, diff < 1 inside the ellipsoid
       diffs = [ abs(self[i]-suite_angles[i])/weights[i] for i in indices ]
       powers = [ d**POWER for d in diffs ]
       return sum(powers) ** (1.0/POWER) # root

    def __str__(self):
        result = "%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f"%\
                 (self.deltam,self.epsilon,self.zeta,self.alpha,self.beta,self.gamma,self.delta)
        return result


class Suite(SuiteAngles):
   """Class that reads dihedral angles out of PDB.Residue objects."""
   def __init__(self):
       SuiteAngles.__init__(self,EMPTY_ANGLES)
       self.resi1 = None
       self.resi2 = None
       self.vecs1 = None
       self.vecs2 = None

   def set_residues(self, resi1, resi2):
       """
       Takes two PDB.Residue objects and calculates suite dihedrals
       from them. If the residues do not apply, a SuiteError is raised.
       """
       rv1 = self.get_vectors(resi1, FIRST_SUITEATOMS)
       rv2 = self.get_vectors(resi2, SECOND_SUITEATOMS)
       self.resi1 = resi1
       self.resi2 = resi2
       self.vecs1 = rv1
       self.vecs2 = rv2
       # calculate torsion angles
       deltam  = dihedral_from_vectors(rv1["C5*"],rv1["C4*"],rv1["C3*"],rv1["O3*"])
       epsilon = dihedral_from_vectors(rv1["C4*"],rv1["C3*"],rv1["O3*"],rv2["P"])
       zeta    = dihedral_from_vectors(rv1["C3*"],rv1["O3*"],rv2["P"],rv2["O5*"])
       alpha   = dihedral_from_vectors(rv1["O3*"],rv2["P"],rv2["O5*"],rv2["C5*"])
       beta    = dihedral_from_vectors(rv2["P"],rv2["O5*"],rv2["C5*"],rv2["C4*"])
       gamma   = dihedral_from_vectors(rv2["O5*"],rv2["C5*"],rv2["C4*"],rv2["C3*"])
       delta   = dihedral_from_vectors(rv2["C5*"],rv2["C4*"],rv2["C3*"],rv2["O3*"])
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
         vec1 = vectors[i]
         vec2 = vectors[i+1]
         dist = vec2-vec1
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

   def confirm_angles(self):
      """checks whether all angles are in the valid range."""
      for a in self.angles:
         if not (0 <= a <= 360):
               raise SuiteAngleError('tangled: invalid angle: %6.3f'%a)

   def confirm(self,check_angles,check_bfactors,check_distances):
       if check_angles: self.confirm_angles()
       if check_bfactors: self.confirm_bfactors()
       if check_distances: self.confirm_distances()
    
   def get_vectors(self,resi, atom_names):
      """Creates a set of Vector objects from a residue object."""
      vectors = {}
      for name in atom_names:
         # new PDB files have ' instead of *
         quotename = name.replace('*',"'")
         if resi.has_id(name):
            vectors[name] = resi[name].get_vector()
         elif resi.has_id(quotename):
            vectors[name] = resi[quotename].get_vector()
         else:
            raise SuiteIncompleteResidueError('atom %s does not exist.'%name)
      return vectors


