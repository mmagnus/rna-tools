#!/usr/bin/env python
#
# StackingCalculator.py
#
# Class that calculates base stacking.
# 
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

"""
A procedure for calculating stacking of RNA nucleotides.

The definition of base stacking from Major & Lemieux 
MC-Annotate paper (JMB 2001, 308, p.919ff):
"Stacking between two nitrogen bases is considered
if the distance between their rings is less
than 5.5 Ang., the angle between the two normals to
the base planes is inferior to 30 deg., and the angle
between the normal of one base plane and the vector
between the center of the rings from the two
bases is less than 40 deg."

There are two classes defined here:
- ResidueVector
- StackingCalculator

The latter class should be used for calculating stacking. There are two 
public methods inside StackingCalculator class that can be used
for calculating stacking:
    
    - process_pdbfile(file_name, chain_id='A') - which runs StackingCalculator
    on the RNA from the 'file_name'.
    The second parameter is optional and has to be set, if the chain ID
    of RNA from PDB file is different than 'A'.
    
"""

import sys
from numpy import array, add, cross, sqrt, arccos
from rna_tools.tools.mini_moderna3.moderna import *
from rna_tools.tools.mini_moderna3.moderna.Constants import NORMAL_SUPPORT, ARCPI

STACKINGS = {
        (True, True): '>>', 
        (True, False): '<<', 
        (False, False): '<>', 
        (False, True): '><', 
        }

# vector placeholder functions
# code snatched from Scientific.Geometry
def angle(vec_a, vec_b):
    cosa = add.reduce(vec_a*vec_b) / \
        sqrt(add.reduce(vec_a*vec_a) * \
        add.reduce(vec_b*vec_b))
    cosa = max(-1., min(1., cosa))
    return arccos(cosa) * ARCPI
    
class StackingInteraction(object):
    """Result from stacking calculation."""
    def __init__(self, resi1, resi2, stack_type):
        """Creates a stacking object."""
        self.resi1 = resi1
        self.resi2 = resi2
        self.type = stack_type
        
    def __repr__(self):
        return "%s %s %s"% \
            (self.resi1.identifier, self.type, self.resi2.identifier)
    
    
class ResidueVector(object):
    """
    Residue class with center vector and normal vector for stacking calculation.
    """
    def __init__(self, residue):
        """
        Creates a dictionary of vectors for each atom from a ModernaResidue.
        """
        self.residue = residue
        self.atoms = {}
        for atom in residue.get_list():
            atom_name = atom.get_fullname().strip().upper()
            self.atoms[atom_name] = residue[atom_name].coord
        self.normal_set = NORMAL_SUPPORT.get(residue.original_base)
        self.normal = None
        self.center = None
            
    def is_valid(self):
        """Checks if all necessary atoms are present."""
        if self.normal_set:
            for name in self.normal_set:
                if name not in self.atoms: 
                    return False
            return True
    
    def calculate_vectors(self):
        """
        Constructs the normal vectors for nucleotide bases.
        Returns a tuple of vectors, the first pointing
        from O to the center of the six-ring of the according base, 
        and the second being the normal 
        vector according to the definition of Major & Thibault 2006.
        Assumes the residue has a complete set of atoms.
        """
        # sum all six atom vectors up to get center point.
        asum = array([0.0, 0.0, 0.0])
        for atomname in self.normal_set:
            asum += self.atoms[atomname]
        self.center = asum / 6.0

        # get two pairs of atoms spanning a plane
        # and calculate the normal vector
        atoma = self.atoms[self.normal_set[1]] - self.atoms[self.normal_set[0]]
        atomb = self.atoms[self.normal_set[3]] - self.atoms[self.normal_set[2]]
        self.normal = cross(atoma, atomb)
        self.normal = self.normal/sqrt(add.reduce(self.normal*self.normal)) 

    def calc_angles(self, rvec):
        """
        Calculates whether the distance and angles between the vectors are OK.
        Returns a tuple of (dist,nn_angle,n1cc_angle,n2cc_angle) or None.
        """
        # calculate the distance between the two ring centers
        ccvec = rvec.center - self.center
        dist = sqrt(add.reduce(ccvec*ccvec)) # vector length
        # check whether the distance is small enough to allow stacking
        if 0.0 < dist < 5.5:
            # check whether the angles are in the allowed range
            nn_angle = angle(self.normal, rvec.normal)
            if (nn_angle < 30 or nn_angle > 150):
                n1cc_angle = angle(self.normal, ccvec)
                n2cc_angle = angle(rvec.normal, ccvec)
                return (dist, nn_angle, n1cc_angle, n2cc_angle)
        return (None, None, None, None)
        
    def get_stacking(self, rvec):
        """
        Returns dictionary with one of the types
        (<<, >>, <>, ><) for the two residues.
        Or None, if they are not stacked.
        """
        distance, nn_ang, n1cc_ang, n2cc_ang = self.calc_angles(rvec)
        if distance and (n1cc_ang < 40 or n1cc_ang > 140 \
                          or n2cc_ang < 40 or n2cc_ang > 140):
            # find out whether the normals are opposed or straight 
            # (pointing in the same direction).
            if nn_ang < 30: 
                straight = True
            elif nn_ang > 150: 
                straight = False
            else: 
                return None # invalid normal angle
            # find out whether base2 is on top of base1
            # calculate whether the normal on base1 brings one closer to base2
            n1c2 = rvec.center - self.center - self.normal
            n1c2dist = sqrt(add.reduce(n1c2*n1c2)) # vector length
            is_up = n1c2dist < distance
            
            stacktype = STACKINGS[(straight, is_up)]
            return StackingInteraction(self.residue, \
                rvec.residue, stacktype)

        
class StackingCalculator:
    """
    Calculates stacking of nucleotide bases according 
    to the definition of Major & Thibault 2006.
        
    Input are residues as parsed by Bio.PDB or Moderna.
    Output are the two residue objects and 
    >> << <> >< stacking codes.
    """
    def get_stacking(self, moderna_struct):
        """
        Loops through all the residues in a ModernaStructure object,
        calls the stacking calculation procedure for all of them.
                
        The method returns list of tuples. Each tuple contains:
        - Residue Object one
        - Residue Object two
        - stacking type of these residues (>>, <<, <> or ><)
        """
        result = []
        rvectors = self.calc_residue_vectors(moderna_struct)
        for record in self.calc_stacking(rvectors):
            if record not in result:
                result.append(record)
        return result
        
    def calc_residue_vectors(self, moderna_struct):
        """
        Precalculates vectors on each residue to make calculations faster.
        """
        rvectors = []
        for residue in moderna_struct:
            rv = ResidueVector(residue) 
            if rv.is_valid():   
                rv.calculate_vectors()  
                rvectors.append(rv)
        return rvectors

    def calc_stacking(self, rvectors):
        """
        Calculates stacking for all residues.
        Generates tuples of (residue1,residue2,stacking_type).
        """
        n_residues = len(rvectors)
        for i in range(n_residues-1):
            resvec1 = rvectors[i]
            for j in range(i+1, n_residues):
                resvec2 = rvectors[j]
                stacking = resvec1.get_stacking(resvec2)
                if stacking:
                    yield stacking

