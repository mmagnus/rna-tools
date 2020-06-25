#!/usr/bin/env python
#
# LIR.py
#
# Calculates LIR values, makes LIRRecord.
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


import math
from math import pi
from numpy import array
from Bio.PDB.vectors import calc_dihedral

# In the production version, the code needs to operate
# on objects rather than files. We need to keep this
# in mind when fixing the API.

#
# CONSTANTS USED FIRST WEEK OF SEP 2009
#
###################### LIR LIBRARY CALCULATION CONSTANTS #######################
### ATOMS FOR CALCULATING VECTORS
#STEM5_AT1 = "C3'"
#STEM5_AT2 = "O3'"
#STEM3_AT1 = "C5'"
#STEM3_AT2 = "C4'"
#
### ATOMS FOR CALCULATING OMEGA DIHEDRALS
#OMEGA5_AT = "N*'"
#OMEGA3_AT = "N*'"
#
### LIR FRAGMENT SUPERPOSITION ATOMS
#STEM5_SUPERPOSITION = ["O3'", "C3'", "C4'"]
#STEM3_SUPERPOSITION = ["C5'", "C4'", "C3'"]
#
### LIR FRAGMENT INSERTION
#STEM5_FRAGMENT = ["C3'", "O3'"]
#STEM5_MODEL = ["C4'", "C5'", "O5'",  "P",  "OP1",  "OP2"]
#STEM5_OTHER = 'model' # indicates whether the rest of atoms (base) should be taken from model or from fragment
#STEM3_FRAGMENT = ["C4'", "C5'", "O5'",  "P",  "OP1",  "OP2"]
#STEM3_MODEL =  ["C3'", "O3'"]
#STEM3_OTHER = 'model' # indicates whether the rest of atoms (base) should be taken from model or from fragment
######################################################################

#
# CONSTANTS TRIED THIRD WEEK OF SEP 2009
#
##################### LIR LIBRARY CALCULATION CONSTANTS #######################
## ATOMS FOR CALCULATING VECTORS
STEM5_AT1 = "C1'"
STEM5_AT2 = "N*"
STEM3_AT1 = "N*"
STEM3_AT2 = "C1'"

## ATOMS FOR CALCULATING OMEGA DIHEDRALS
OMEGA5_AT = "C4'"
OMEGA3_AT = "C4'"
#####################################################################


class LirRecord(object):
    """
Represents a single fragment as it is stored in the database

fr_length    - integer
structure      - PDB ID
chain          - chain ID (char)
preceding_resi - residue number immediately before the fragment
sequence       -  Structure object

    """
    def __init__(self, fr_length=None, structure=None, chain=None, preceding_resi=None, following_resi=None, \
                 sequence=None, sequence_anchor=None, secstruc=None, x=None, y=None, dist_anchor=None, beta=None, gamma=None, omega5=None, omega3=None, \
                 P_dist=None, O5p_dist=None, C5p_dist=None, C4p_dist=None, C3p_dist=None, O3p_dist=None, O2p_dist=None, C1p_dist=None, N_dist=None):
        #TODO: dictionary
        self.fr_length = fr_length
        self.structure = structure
        self.chain = chain
        self.preceding_residue = preceding_resi
        self.following_residue = following_resi
        self.sequence = sequence
        self.sequence_anchor = sequence_anchor
        self.secstruc = secstruc
        self.x = x
        self.y = y
        self.dist_anchor = dist_anchor
        self.beta = beta
        self.gamma = gamma
        self.omega5 = omega5
        self.omega3 = omega3
        
        self.distances = array([P_dist, O5p_dist, C5p_dist, C4p_dist, C3p_dist, O3p_dist, O2p_dist, C1p_dist, N_dist])

    def __nonzero__(self):
        """A LirRecord object is always not empty"""
        return True

    def __str__(self):
        """
        """
        result = '-'*50
        result += '\nFragment length: %i'%self.fr_length
        result += 'Structure: %s '%(self.structure)+'  Chain: %s'%self.chain+'  \nSequence: %s' %self.sequence.seq_with_modifications
        result += '\nPreceding residue: %s'%self.preceding_residue+'Following residue: %s'%self.following_residue
        # result += '\nx: %5.2f'%self.x+ '  y: %5.2f'%self.y+ '  Dist anchor: %5.2f'%self.dist_anchor+ '  Beta: %5.2f'%self.beta+'  Dihedral: %5.2f'%self.gamma
        result += '-'*50
        return result


    def get_list(self):
        """Returns a list which contains all fields of a LirRecord"""
        result = [
            self.fr_length, self.structure, self.chain,
            self.preceding_residue, self.following_residue,  self.sequence, self.sequence_anchor, self.secstruc]
        result += [d for d in self.distances]        
        return result


    def get_txt(self, separator):
        """
        Returns a string containing all fields of a LirRecord
        The separator character between values can be defined 
        by the user (e.g. '|')
        """
        seq = self.sequence.seq_with_modifications
        seq_anchor = self.sequence_anchor.seq_with_modifications
        elem_list = [str(elem) for elem in self.get_list()]
        elem_list[5] = seq
        elem_list[6] = seq_anchor
        result = separator.join(elem_list)
        return result


class Lir:
    """
Calculates parameters of single fragment.
These are dihedral angles (gamma, omega5, omega3),
flat angle beta, x and y values, anchor distance.
See : LIP - Loops In Proteins. Michalska et al.
    """
    def __init__(self, anchor5,anchor3):
        """
        """
        self.anchor5=anchor5
        self.anchor3=anchor3
        self.set_vectors()
        self.calculate_atom_distances()
        self.calculate_lir_values()
        
    def set_vectors(self):
        """
        Getting vectors for 2 main atoms from anchor residues.
        # Stem5 ---> P, C4';    Stem3 ---> C4', O3';
        """
        # What can go wrong: such atoms may not be present in residue.
        self.vP5=self.anchor5.get_atom_vector(STEM5_AT1)
        self.vC5=self.anchor5.get_atom_vector(STEM5_AT2)
        self.vC3=self.anchor3.get_atom_vector(STEM3_AT1)
        self.vO3=self.anchor3.get_atom_vector(STEM3_AT2)
        # Getting vectors for base atom from anchor residues.
        self.vN5 = self.anchor5.get_atom_vector(OMEGA5_AT)
        self.vN3 = self.anchor3.get_atom_vector(OMEGA5_AT)
            
        self.vP5C5=self.vC5-self.vP5
        self.vP5C3=self.vC3-self.vP5

        vnP5C5=self.vP5C5/self.vP5C5.norm() # length
        vnP5C3=self.vP5C3/self.vP5C3.norm() # length
        self.cosC5P5C3=vnP5C5*vnP5C3


    def calculate_dihedrals(self):
        """
        Calculates the dihedral angles: 
        """
        # create positive dihedrals [0..2pi]
        posi = lambda x:x<0 and 2*pi+x or x
        self.gamma  = posi(calc_dihedral(self.vP5,self.vC5,self.vC3,self.vO3))
        self.omega5 = posi(calc_dihedral(self.vO3,self.vC3,self.vC5,self.vN5))
        self.omega3 = posi(calc_dihedral(self.vP5,self.vC5,self.vC3,self.vN3))
    

    def calculate_beta(self):
        """
        Calculate flat angle: C4' (anchor5) --- C4' (anchor3) --- O3' (anchor3)
        """
        vC5C3 = self.vC3-self.vC5
        vC3O3 = self.vC3-self.vO3
        self.beta = vC5C3.angle(vC3O3)

    def calculate_x_value(self):
        """
        Returns the length of the distance vector between both anchors
        parallel to the first anchor.
        """
        self.x=abs((self.cosC5P5C3*self.vP5C3.norm())-self.vP5C5.norm()) # norm() is length()

    def calculate_y_value(self):
        """
        Returns the length of the distance vector between both anchors
        perpendicular to the first anchor.
        """
        self.y=abs(math.sin(math.acos(self.cosC5P5C3))*self.vP5C3.norm()) # norm() is length()

    def calculate_dist_anchor(self):
        """
        Calculates the Euclidean distance of the two anchor residues.
        (the start and end points of the fragment).
        """
        self.dist_anchor = math.sqrt(self.x**2+self.y**2)
    
    def calculate_atom_distances(self):
        self.P_dist = self.anchor5['P'] - self.anchor3['P']
        self.O5p_dist = self.anchor5["O5'"] - self.anchor3["O5'"]
        self.C5p_dist = self.anchor5["C5'"] - self.anchor3["C5'"]
        self.C4p_dist = self.anchor5["C4'"] - self.anchor3["C4'"]
        self.C3p_dist = self.anchor5["C3'"] - self.anchor3["C3'"]
        self.O3p_dist = self.anchor5["O3'"] - self.anchor3["O3'"]
        self.O2p_dist = self.anchor5["O2'"] - self.anchor3["O2'"]
        self.C1p_dist = self.anchor5["C1'"] - self.anchor3["C1'"]
        self.N_dist = self.anchor5['N*'] - self.anchor3['N*']
        
    def calculate_lir_values(self):
        self.calculate_dihedrals()
        self.calculate_beta()
        self.calculate_x_value()
        self.calculate_y_value()
        self.calculate_dist_anchor()
