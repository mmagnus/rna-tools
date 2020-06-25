#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# PuckerCalculator.py
#
# Calculates the ribose pucker in nucleic acid structures.
#
__author__ = "Tomasz Osinski"
__copyright__ = "Genesilico 2008"
__credits__ = ["Kristian Rother", "Raphael Bauer", "Marcin Domagalski", \
    "Magdalena Rother", "Janusz Bujnicki", "Marie Curie"]
__license__ = "GPL"
__status__ = "Production"

from Bio.PDB.Vector import calc_dihedral        
from math import pi, sin, atan, degrees

FURANOSEATOMS = ["O4'", "C1'", "C2'", "C3'", "C4'"]

THETA_RANGE = (
    ( 72,108, "C3'-endo"),
    (108,144, "C4'-exo"),
    (144,180, "O4'-endo"),
    (180,216, "C1'-exo"),
    (216,252, "C2'-endo"),
    (252,298, "C3'-exo"),
    (298,324, "C4'-endo"),
    (324,360, "O4'-exo"),
    (  0, 36, "C1'-endo"),
    ( 36, 72, "C2'-exo"),
    )

class PuckerError(Exception): pass

def check_positive(angle):
    """Turns angle positive"""
    if angle < 0: 
        angle += 360.0
    return angle

class PuckerCalculator:
    """Class that parses ribose puckers out of PDB.Residue objects."""
    def __init__(self):
        """Initializes the calculator."""
        self.resi = None
        self.vectors = {}
        self.angles = [0.0] * 5
        self.pucker_angle = 0.0

    def __repr__(self):
        """Returns recently calculated dihedrals."""
        return "%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f"% \
                 (self.angles[0], self.angles[1], self.angles[2], \
                    self.angles[3], self.angles[4], self.pucker_angle)

    def _calc_vectors(self):
        """Creates a set of Vector objects from a residue object."""
        self.vectors = {}
        for name in FURANOSEATOMS:
            self.vectors[name] = self.resi[name].get_vector()

    def _dihedral(self, vec1, vec2, vec3, vec4):
        """Calculates torsion and makes sure it is between -180 and +180"""
        torsion = calc_dihedral(vec1, vec2, vec3, vec4)
        if torsion > 180:
            return -(360-torsion)
        else:
            return torsion
        
    def _calc_torsions(self):
        """Calculate torsion angles"""
        rv1 = self.vectors
        self.angles = (
            self._dihedral(rv1["C4'"], rv1["O4'"], rv1["C1'"], rv1["C2'"]),
            self._dihedral(rv1["O4'"], rv1["C1'"], rv1["C2'"], rv1["C3'"]),
            self._dihedral(rv1["C1'"], rv1["C2'"], rv1["C3'"], rv1["C4'"]),
            self._dihedral(rv1["C2'"], rv1["C3'"], rv1["C4'"], rv1["O4'"]),
            self._dihedral(rv1["C3'"], rv1["C4'"], rv1["O4'"], rv1["C1'"]),
            )

    def _calc_pucker_angle(self):
        """Determines the pucker angle."""
        self.pucker_angle = degrees(atan(((self.angles[2]+self.angles[4]) \
                        - (self.angles[1]+self.angles[3])) \
                        /(2*self.angles[0]*(sin(36*pi/180.0)\
                        +sin(72*pi/180.0)))))
        if self.angles[0] < 0:
            self.pucker_angle += 180
        if self.pucker_angle < 0: 
            self.pucker_angle += 360
        #self.angles = map(check_positive, self.angles)
    
    def find_pucker(self):
        """Determines the pucker type."""
        for min_ft, max_ft, pucker in THETA_RANGE:
            if min_ft < self.pucker_angle <= max_ft:
                return pucker
        
    def get_pucker(self, resi):
        """Returns the pucker for a PDB.Residue object."""
        self.resi = resi
        self._calc_vectors()
        self._calc_torsions()
        self._calc_pucker_angle()
        return self.find_pucker()




