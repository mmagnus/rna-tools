#!/usr/bin/env python
#
# FCCDLoopCloser.py
#
# For backbone fixing.
#
# http://iimcb.genesilico.pl/moderna/
#

"""
FCCD Loop Closing Algorithm

concept by M.Boomsma and T.Hamelryck
implemented by Kristian Rother and Magdalena Rother.

when referring to this algorithm, please refer to

Wouter Boomsma and Thomas Hamelryck,
Full cyclic coordinate descent: solving the protein loop closure problem
in Calpha space,
BMC Bioinformatics 2005, 6:159doi:10.1186/1471-2105-6-159.
"""

__author__ = "Kristian Rother, Magdalena Rother, Tomasz Puton"
__copyright__ = "Copyright 2008, The Moderna Project"
__license__ = "GPL"
__credits__ = ["Janusz Bujnicki"]
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"


from math import sqrt
from numpy import array, dot, zeros, transpose, eye
from numpy.linalg import linalg

S = eye(3)
S[2, 2] = -1

class FCCDError(Exception): pass

class FCCDLoopCloser:
    """
    Class for running the FCCD loop closure algorithm
    on a list of fixed and a set of moved atoms.
    """
    def __init__(self, moving, fixed):
        """
        moving - list of 6+ atoms that are moved.
        fixed - list of three atoms that stay the same.
        All atoms should be Bio.PDB.Atom objects.
        """
        if len(moving) < 6:
            raise FCCDError("""Moving should have at least length 6
    (3 atoms at start, plus 3 to match with fixed.""")
        if len(fixed) != 3:
            raise FCCDError(""""Fixed should have length 3
    (atoms at the end to be closed)""")
        self.moving_atoms = moving
        self.moving = [m.get_vector() for m in moving]
        self.fixed = [f.get_vector() for f in fixed]
        # Coordinates along COLUMNS
        self.moving_coords = zeros((3, 3), 'd')
        self.fixed_coords = zeros((3, 3), 'd')
        
    def calc_rmsd(self):
        """Returns RMSD fit of last 3 moving vectors to fixed vectors."""
        rmsd = 0.0
        for i_vec in range(1, 4):
            dist = self.moving[-i_vec] - self.fixed[-i_vec]
            rmsd += dist.norm()**2
        return sqrt(rmsd/3.0)
        
    def copy_vectors_to_atoms(self):
        """Copies the current coordinates to atom objects."""
        for vec, atom in zip(self.moving, self.moving_atoms[:-3]):
            atom.coord = array([vec[0], vec[1], vec[2]])
        
    def get_moving_coords(self, center):
        """
        move to pivot origin
        Prepares arrays with the shifted coordinates
        of the three last atoms of both chains. 
        """
        for i_vec in range(0, 3):
            index = -(3 - i_vec)
            vec = self.moving[index] - center
            self.moving_coords[:, i_vec] = vec.get_array()
        return self.moving_coords

    def get_fixed_coords(self, center):
        """
        move to pivot origin
        Prepares arrays with the shifted coordinates
        of the three last atoms of both chains. 
        """
        for i_vec in range(0, 3):
            vec = self.fixed[i_vec] - center
            self.fixed_coords[:, i_vec] = vec.get_array()
        return self.fixed_coords

    def accept_rotation(self):
        """
        Abstract method that could be used to apply
        angle constraints or Metropolis criterion.
        """
        return True

    def apply_rotation(self, rot_matrix, i_vec, center):
        """Adjusts the coordinates"""
        for j_vec in range(i_vec+1, len(self.moving)):
            vec = self.moving[j_vec] - center
            vec = vec.left_multiply(rot_matrix)
            vec = vec + center
            self.moving[j_vec] = vec


    def optimize_vector(self, i_vec):
        """
        Optimizes torsion and bond angles by rotation around
        one atom i.
        1. moves coordinate origin to that atom.
        2. generate a rotation matrix from a singular value decompositon
        3. rotate all atoms after i.
        """
        center = self.moving[i_vec]
        moving_coords = self.get_moving_coords(center)
        fixed_coords = self.get_fixed_coords(center)
        # Do singular value decomposition
        a = dot(fixed_coords, transpose(moving_coords))
        u, d, vt = linalg.svd(a)
        # Check reflection
        if (linalg.det(u) * linalg.det(vt))<0:
            u = dot(u, S)
        # Calculate rotation
        rot_matrix = dot(u, vt)
        # Apply rotation
        if self.accept_rotation():
            self.apply_rotation(rot_matrix, i_vec, center)


    def run_fccd(self, threshold=0.2, maxit=100):
        """
        The moving vectors are changed until its last three elements
        overlap with the fixed ones with a RMSD smaller than threshold.
        """
        n_it = 0
        while n_it < maxit:
            for i_vec in range(2, len(self.moving) - 2):
                self.optimize_vector(i_vec)
                rmsd = self.calc_rmsd()
                if rmsd < threshold:
                    return "RMSD threshold reached", rmsd, n_it
                n_it += 1
        return "Stop - max iterations reached", rmsd, n_it
            



    

    



