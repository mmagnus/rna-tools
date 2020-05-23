"""
implementation of Moderna algorithm for detecting stackings
"""
import simplejson as json
from xml.dom import minidom

import sys, os, re

import numpy as np
import itertools
import string
import gzip
__email__ = "gchojnowski@genesilico.pl"

from operator import itemgetter
from string import strip

from cStringIO import StringIO

from Bio import PDB
from scipy.spatial import KDTree
from optparse import OptionParser


# *************************************************************************************************
# *************************************************************************************************
# *************************************************************************************************
# *************************************************************************************************
# Moderna stacking detector


#__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
#__copyright__ = "Copyright 2008, The Moderna Project"
#__credits__ = ["Janusz Bujnicki"]
#__license__ = "GPL"
#__version__ = "1.5.0"
#__maintainer__ = "Magdalena Rother"
#__email__ = "mmusiel@genesilico.pl"
#__status__ = "Production"
# ADAPTED for CCTBX objects from: StackingCalculator.py
# 27.06.2011 gchojnowski@genesilico.pl

#from numpy import array, add, cross, sqrt, arccos
#import math

# between purines and pyrimidines, the normals are reversed, because
# the rings are reversed with respect to the helix axis.
NORMAL_SUPPORT = {
        'C':['N1','C2','N3','C4','C5','C6'],
        'U':['N1','C2','N3','C4','C5','C6'],
        'T':['N1','C2','N3','C4','C5','C6'],
        'G':['N1','C2','C4','N3','C5','C6'],
        'A':['N1','C2','C4','N3','C5','C6'],
        }


STACKINGS = {
        (True, True): '>>',
        (True, False): '<<',
        (False, False): '<>',
        (False, True): '><',
        }

ARCPI = 180.0/np.pi



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
# *************************************************************************************************


class ResidueVector:

    def __init__(self, residue):
        """Creates a dictionary of vectors for each atom from a ModernaResidue."""
        self.residue = residue
        self.atoms = {}

        for atom in residue:
            atom_name = atom.name.strip().upper()
            self.atoms[atom_name] = np.array(atom.get_coord())

        self.resn = residue.resname.strip()
        self.normal_set = NORMAL_SUPPORT.get(self.resn.upper())
        self.normal = None
        self.center = None

    # vector placeholder functions
    # code snatched from Scientific.Geometry
    def angle(self, vec_a, vec_b):
        cosa = np.add.reduce(vec_a*vec_b) / \
            np.sqrt(np.add.reduce(vec_a*vec_a) * \
            np.add.reduce(vec_b*vec_b))
        cosa = max(-1., min(1., cosa))
        return np.arccos(cosa) * ARCPI


    def is_valid(self):
        """Checks if all necessary atoms are present."""
        if self.normal_set:
            for name in self.normal_set:
                if not self.atoms.has_key(name):
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
        asum = np.array([0.0, 0.0, 0.0])
        for atomname in self.normal_set:
            if not self.atoms.has_key(atomname):
                self.normal = None
                return
            asum += self.atoms[atomname]
        self.center = asum / 6.0

        # get two pairs of atoms spanning a plane
        # and calculate the normal vector
        atoma = self.atoms[self.normal_set[1]] - self.atoms[self.normal_set[0]]
        atomb = self.atoms[self.normal_set[3]] - self.atoms[self.normal_set[2]]
        self.normal = np.cross(atoma, atomb)
        self.normal = self.normal/np.sqrt(np.add.reduce(self.normal*self.normal))



    def calc_angles(self, rvec):
        """
        Calculates whether the distance and angles between the vectors are OK.
        Returns a tuple of (dist,nn_angle,n1cc_angle,n2cc_angle) or None.
        """
        # calculate the distance between the two ring centers
        ccvec = rvec.center - self.center
        dist = np.sqrt(np.add.reduce(ccvec*ccvec)) # vector length
        # check whether the distance is small enough to allow stacking
        if 0.0 < dist < 5.5:
            # check whether the angles are in the allowed range
            nn_angle = self.angle(self.normal, rvec.normal)
            if (nn_angle < 30 or nn_angle > 150):
                n1cc_angle = self.angle(self.normal, ccvec)
                n2cc_angle = self.angle(rvec.normal, ccvec)
                return (dist, nn_angle, n1cc_angle, n2cc_angle)
        return (None, None, None, None)


    def get_stacking(self, rvec):
        """
        Returns dictionary with one of the types
        (<<, >>, <>, ><) for the two residues.
        Or None, if they are not stacked.
        """
        if self.normal is None or rvec.normal is None:
            return None
        distance, nn_ang, n1cc_ang, n2cc_ang = self.calc_angles(rvec)
        #print n1cc_ang, n2cc_ang
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
            n1c2dist = np.sqrt(np.add.reduce(n1c2*n1c2)) # vector length
            is_up = n1c2dist < distance

            stacktype = STACKINGS[(straight, is_up)]
            return stacktype#self.residue, rvec.residue, stacktype




# *************************************************************************************************
# *************************************************************************************************
# *************************************************************************************************
# *************************************************************************************************





class contacts:

    def __init__(self, structure):
        self.structure = structure

    # -------------------------------------------------------------------------

    def calc_stacking(self, resi_1, resi_2):
        residue_vector1 = ResidueVector(resi_1)
        residue_vector2 = ResidueVector(resi_2)

        if not residue_vector1.normal_set or not residue_vector2.normal_set:
            return None
        else:
            residue_vector1.calculate_vectors()
            residue_vector2.calculate_vectors()
            return residue_vector1.get_stacking(residue_vector2)

    # -------------------------------------------------------------------------

    def run(self,use_fr3d_format=False):
        use_fr3d_format_old = False
        if use_fr3d_format_old:
            print "# stacking"
        elif use_fr3d_format:
            pass
        else:
            print "Stackings -------------------------------------------------------"
        residues = [r for r in self.structure.get_residues()]
        n = len(residues)
        r_coord = []
        for r in residues:
            p = None
            for a in r:
                if p is None:
                    p = a.get_coord()
                if a.name.strip() == "P":
                    p = a.get_coord()
            r_coord.append(p)

        tree = KDTree(r_coord)
        residues_name = [ "%s%s" % (r.get_parent().get_id(), r.get_id()[1]) for r in residues]
        n_types = [r.get_resname().strip() for r in residues]
        num = 0
        for i in xrange(n):
            for j in tree.query(r_coord[i], k=100)[1]:
                if i>=j or j>=len(residues):
                    continue
                stacking = self.calc_stacking(residues[i], residues[j])
                if stacking is not None:
                    if use_fr3d_format:
                        # example:
                        # "0_A_1_U","s35","0_A_2_G"
                        num += 1
                        n_type_i = n_types[i]
                        n_type_j = n_types[j]
                        chain_i = residues_name[i][0]
                        num_i = residues_name[i][1:]
                        chain_j = residues_name[j][0]
                        num_j = residues_name[j][1:]
                        print '"1_%s_%s_%s","%s","1_%s_%s_%s"' % (chain_i,num_i,n_type_i,stacking,chain_j,num_j,n_type_j)
                    elif use_fr3d_format_old:
                        # example:
                        # 2210 G1003A(A) - G1003(A) -    s53  -    0
                        num += 1
                        n_type_i = n_types[i]
                        n_type_j = n_types[j]
                        print "%-4d %s(%s) - %s(%s) - %s - 0" % (num, residues_name[i],n_type_i,residues_name[j],n_type_j,stacking)
                    else:
                        print "%s-%s : %s" % (residues_name[i],residues_name[j],stacking)

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""some experiments with classifier generation""")
    parser.add_option("--use-fr3d-format", dest="use_fr3d_format", action='store_true',
                  help="use fr3d output format",default=False)
    parser.add_option("--use-fr3d-format-old", dest="use_fr3d_format_old", action='store_true',
                  help="use fr3d output format (OLD)",default=False)

    (options, args)  = parser.parse_args()
    return (parser, options, args)


def process_structures(filename,use_fr3d_format):
    if re.match("^.*.gz$",filename):
        f = gzip.open(filename)
    else:
        f = filename
    parser = PDB.PDBParser()
    structure = parser.get_structure("c", f)
    contacts_obj = contacts(structure)
    contacts_obj.run(use_fr3d_format)

if __name__=="__main__":
    (_parser,options,args) = parse_args()
    process_structures(args[0],options.use_fr3d_format)




