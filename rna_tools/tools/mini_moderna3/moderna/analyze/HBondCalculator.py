#!/usr/bin/env python
#
# HBondCalculator.py
#
# Calculates hydrogen bonds between residues.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@genesilico.pl"
__status__ = "Prototype"
# -*- coding: cp1250 -*-

from Bio.PDB.vectors import calc_angle
from math import degrees
from Bio.PDB.vectors import Vector
from rna_tools.tools.mini_moderna3.moderna.Constants import MAX_C1_DIST, MAXD1, MAXD2
from rna_tools.tools.mini_moderna3.moderna.Constants import BACKBONE_RIBOSE_ATOMS
#TODO: bifurcated bonds
#TODO: different parameter sets
#TODO: h-bond tolerance depending on b-factor

def sq3_norm(a):
    return a.norm()
    
class HBond:
    def __init__(self,donor,acceptor,dist):
        self.donor = donor
        self.acceptor = acceptor
        self.dist = dist
        self.donor_resi = donor.parent
        self.acceptor_resi = acceptor.parent

    def is_primary(self):
        """Returns True if hbond fulfills criteria for primary hbond:
It should be a H-bond that occurs between donor and acceptor
atoms, both of which belong to a base, and with the distance
between the two atoms (N or O) <3.4A.
"""
        if self.donor.name in BACKBONE_RIBOSE_ATOMS: return False
        if self.acceptor.name in BACKBONE_RIBOSE_ATOMS: return False
        if self.donor.name[0] == 'C' or self.acceptor.name[0] == 'C': return False
        if self.dist >= 3.4: return False
        return True

    def is_secondary(self):
        """Returns True if hbond fulfills criteria for secondary hbond:
The second H-bond is determined with extended distances.
The distance between donor and acceptor
for base (N, O) - base (N, O)<3.75A,
for base (N, O) - base (CH) <3.9A
for base (N, O) - ribose (O20) <3.75A
"""
        cutoff = 3.75
        if self.donor.name[0] == 'C': cutoff = 3.9
        if self.dist >= cutoff: return False
        return True
        
    def is_tertiary(self):
        """Returns True if hbond fulfills criteria for tertiary interaction hbond."""
        return True

    def __repr__(self):
        return "%s-%s"%(self.donor.name,self.acceptor.name)
        #return "%s%s-%s%s"%(self.donor.parent.resname,self.donor.name,self.acceptor.parent.resname,self.acceptor.name)


class HBondParameterCalculator:
    """Class that performs the vector calculations for an hbond."""
    def __init__(self, donor,  hydrogen, acceptor, acc_support):
        """Takes four sets of coord vectors."""
        self.donor = Vector(donor)
        self.acceptor = Vector(acceptor)
        self.hydrogen = Vector(hydrogen)
        self.acc_support = Vector(acc_support)
        # helper vectors
        self.dh = self.hydrogen-self.donor
        self.ha = self.hydrogen-self.acceptor
        self.acs = self.acc_support-self.acceptor
        # calculate distances and angles
        
    def _get_d1(self):
        conn = self.acceptor-self.hydrogen
        return conn.norm()

    def _get_d2(self):
        conn = self.donor-self.acceptor
        return conn.norm()
        
    def _get_alpha(self):
        # SLOWER
        #return PDB.calc_angle(Vector(self.donor), Vector(self.hydrogen), Vector(self.acceptor))* 180.0/pi
        return degrees(calc_angle(self.donor, self.hydrogen, self.acceptor))
        
    def _get_beta(self):
        return degrees(calc_angle(self.hydrogen, self.acceptor, self.acc_support))
    
    def _get_gamma(self):
        return degrees(calc_angle(self.donor, self.hydrogen, self.acc_support ))
    
    d1 = property(_get_d1)
    d2 = property(_get_d2)
    alpha = property(_get_alpha)
    beta = property(_get_beta)
    gamma = property(_get_gamma)
    
    def is_valid(self):
        """
        Returns boolean if the parameters are within range        
            Criteria set I (default values) according to Suehnel/HBExplore
            d1 < 3.9 A, d2 < 2.5 A , a > 90, b > 90, c > 90. 
            There is a mistake in the figure. If you read the text you find 
            that d1 and d2 are mixed up! It should be d1<2.5 and d2<3.9
        """
        if self.d2<MAXD2 and self.d1 < MAXD1 \
            and self.alpha > 90 and self.beta>90 and self.gamma>90:
                return True
                
    def __repr__(self):
        return """alpha:%5.1f  beta:%5.1f  gamma:%5.1f  d1:%5.2f  d2:%5.2f"""%(self.alpha, self.beta, self.gamma, self.d1, self.d2)


class HBondCalculator:
    """
    Calculates hydrogen bonds out of pairs of Residue objects.
    (moderna.ModernaResidue or Bio.PDB.Residue)
    """
    def calc_hbond_list(self,res1,res2):
        """
        Calculates and returns a list of HBond objects.
        res1, res2 - Residue objects        
        """
        # shortcut: compare c1'-c1' distance first
        c1dist = res1["C1'"]-res2["C1'"]
        if c1dist > MAX_C1_DIST: return []
        
        hbonds = []
        for acceptor in res1.get_hbond_acceptors():
            for donor in res2.get_hbond_donors():
                hb = self.calc_hbond(donor,acceptor)
                if hb: hbonds.append(hb)
        for acceptor in res2.get_hbond_acceptors():
            for donor in res1.get_hbond_donors():
                hb = self.calc_hbond(donor,acceptor)
                if hb: hbonds.append(hb)
        return hbonds
        

    def get_acceptor_support(self, acceptor):
        """
        Returns coordinates of the support point behind the acceptor atom.
        According to Fig.1 in the Suehnel HBExplore pare, this is A1 if there is
        only one neighbor, and Am if there are two.
        Atoms with 3 or more neighbors will be rejected.
        """
        acc_support = None
        neighbors = acceptor.parent.get_neighbors(acceptor)
        if len(neighbors)==1:
            acc_support = neighbors[0].coord
        elif len(neighbors) == 2:
            acc_support = (neighbors[0].coord + neighbors[1].coord ) * 0.5
        return acc_support
                
    def calc_hbond(self,donor,acceptor):
        """
        Returns a HBond object for two atoms, or None.
        donor, acceptor - Bio.PDB.Atom objects
        """
        dist = donor-acceptor
        if dist < MAXD2: #TODO: refactor this one out (checked twice)
            acc_support = self.get_acceptor_support(acceptor)
            if acc_support != None:
                for hydrogen in donor.parent.get_donor_hydrogens(donor):
                    params = HBondParameterCalculator(donor.coord, hydrogen, acceptor.coord, acc_support)
                    if params.is_valid():
                        return HBond(donor,acceptor,dist)
