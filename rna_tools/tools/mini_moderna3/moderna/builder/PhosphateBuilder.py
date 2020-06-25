#!/usr/bin/env python
#
# PhosphateBuilder.py
#
# Module for reconstructing P, O5', OP1, OP2 atoms
# e.g. after loop insertions.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Prototype"

"""
Constructs backbone atoms that satisfy standard distance, angle, and torsion 
constraints. P and O5' atoms are constructed given that all other backbone atoms are present.

The procedure first places the two atoms at 1/3 (P) and 1/3 (O5') of the O3'-C5' vector.
Then the coordinates are minimized, using random displacements iteratively,
and scoring distance and angle deviations.

The minimization algorithm is very crude and not 100% reliable, but it works a lot
better than the previous approach (pure fragment replacement).

"""

from Bio.PDB import Atom
from numpy import array, ndarray
from random import random
from rna_tools.tools.mini_moderna3.moderna.analyze.GeometryParameters import GeometryStandards
from rna_tools.tools.mini_moderna3.moderna.RNAChain import RNAChain
from rna_tools.tools.mini_moderna3.moderna.ModernaSuperimposer import ModernaSuperimposer
from Bio.PDB.vectors import Vector, calc_angle, calc_dihedral
from rna_tools.tools.mini_moderna3.moderna.Constants import DATA_PATH
from rna_tools.tools.mini_moderna3.moderna.builder.CoordBuilder import build_coord, get_ref_matrix
from rna_tools.tools.mini_moderna3.moderna.analyze.RNASuites import BETA, DEFAULT_TORSIONS
from Bio.PDB import NeighborSearch
import math

gs = GeometryStandards()
O3_P_DIST  = gs.get_standard("X:O3',X+1:P")
P_O5_DIST  = gs.get_standard("X:P,X:O5'")
O5_C5_DIST = gs.get_standard("X:O5',X:C5'")
A1 = gs.get_standard("X:C3',X:O3',X+1:P")
A2 = gs.get_standard("X:O3',X+1:P,X+1:O5'")
A3 = gs.get_standard("X:P,X:O5',X:C5'")
A4 = gs.get_standard("X:O5',X:C5',X:C4'")

P_OP_DIST = 1.485
P_OP_ANGLE = 109.0
OP1_TORSION = 54.1
OP2_TORSION = -178.0
OP3_TORSION = -62.1

STD_ANGLES = array([A2, A3, A4])
ANGLE_WEIGHT = 1.0/180.0
MAKE_DEGREES = 180/math.pi

class PhosphateBuilder(object):
    """
    uses a simplistic approach to construct missing
    C5' and P atoms.
    """
    def __init__(self,resi1,resi2):
        self.r1 = resi1
        self.r2 = resi2
        self.p = None
        self.o5 = None
        self.o3 = self.r1["O3'"]
        self.o3v = self.o3.get_vector()
        self.c3 = self.r1["C3'"]
        self.c3v = self.c3.get_vector()
        self.c4 = self.r1["C4'"]
        self.c4v = self.c4.get_vector()
        self.c5 = self.r2["C5'"]
        self.c5v = self.c5.get_vector()
        self.c42 = self.r2["C4'"]
        self.c42v = self.c42.get_vector()
    
    def build(self):
        """Runs the full phosphate construction procedure."""
        self.remove_old_linker()
        self.construct_initial_atoms()
        self.minimize_atoms()
        self.add_op12()


    def remove_old_linker(self):
        """Removes eventually existing backbone atoms."""
        if self.r2.child_dict.has_key("P"):
            self.r2.detach_child("P")
        if self.r2.child_dict.has_key("O5'"):
            self.r2.detach_child("O5'")
        self.remove_op1_op2()
        
    def remove_op1_op2(self):
        if self.r2.child_dict.has_key("OP1"):
            self.r2.detach_child("OP1")
        if self.r2.child_dict.has_key("OP2"):
            self.r2.detach_child("OP2")
        if self.r2.child_dict.has_key("OP3"):
            self.r2.detach_child("OP3")

    def construct_initial_atoms(self):
        """Finds positions in between the O3' and C5'"""
        vec = self.c5.coord - self.o3.coord
        coord_p = self.o3.coord + vec*0.333333
        coord_o5 = self.o3.coord + vec*0.666666
        self.p = Atom.Atom("P", coord_p, 0.0,1.0,' '," P",1, element="P")
        self.o5 = Atom.Atom("O5'", coord_o5, 0.0,1.0,' '," O5'",1, element="O")
        self.r2.add(self.p)
        self.r2.add(self.o5)
        
    def __str__(self):
        d1 = abs(O3_P_DIST - float(self.o3-self.p))
        d2 = abs(P_O5_DIST - float(self.p-self.o5))
        d3 = abs(O5_C5_DIST - float(self.o5-self.c5))
        a1 = math.degrees(calc_angle(self.c3v,self.o3v,self.p.get_vector()))
        a2 = math.degrees(calc_angle(self.o3v,self.p.get_vector(),self.o5.get_vector()))
        a3 = math.degrees(calc_angle(self.p.get_vector(),self.o5.get_vector(),self.c5.get_vector()))
        a4 = math.degrees(calc_angle(self.o5.get_vector(),self.c5v,self.c42v))
        a1 = abs(a1-A1)
        a2 = abs(a2-A2)
        a3 = abs(a3-A3)
        a4 = abs(a4-A4)
        result = "%5.2f\t%5.2f\t%5.2f\t"%(d1, d2, d3)
        result += "%5.2f\t%5.2f\t%5.2f\t%5.2f\t"%(a1, a2, a3, a4)
        result += "total %5.2f\n"%self.get_score()
        return result
        

    def get_score(self, max_score=0.0):
        """Returns the square deviation from reference values."""
        d2 = P_O5_DIST - float(self.p-self.o5)
        d3 = O5_C5_DIST - float(self.o5-self.c5)
        dscore = d2*d2 + d3*d3
        if dscore > max_score: return dscore
        #a2 = calc_angle(self.o3v,self.pv,self.o5v)* MAKE_DEGREES
        a3 = calc_angle(self.pv,self.o5v,self.c5v)* MAKE_DEGREES
        a4 = calc_angle(self.o5v,self.c5v,self.c42v)* MAKE_DEGREES
        # calc difference to standard values
        #a2 = (a2-A2)*ANGLE_WEIGHT
        a3 = (a3-A3)*ANGLE_WEIGHT
        a4 = (a4-A4)*ANGLE_WEIGHT
        # using a small numpy array is slower
        #ang = array((a2, a3, a4))
        #ang = ang * MAKE_DEGREES
        #ascore = ang-STD_ANGLES
        #ascore *= ANGLE_WEIGHT
        #ascore = sum(ascore * ascore)
        ascore = a3 *a3 + a4 * a4 # a2*a2 + 
        return dscore + ascore


    def get_random_torsion(self,center,step):
        """returns a random dihedral close to a reference value."""
        half = step * 0.5
        torsion = center + random()*step - half
        while torsion > 180: torsion -= 360
        while torsion < -180: torsion += 360
        return torsion
        
        
    def minimize_atoms(self):
        """Performs crude random minimization of dihedral angles."""
        epsilon = 0.02
        i = 0
        maxcyc = 50
        step = 120
        delta = 100 + epsilon
        MAX_J = 5
        MAX_K = 5
        tor1, tor2 = 0.0, 0.0
        matrix1 = get_ref_matrix(self.c4v, self.c3v, self.o3v)
        while delta > epsilon and i < maxcyc:
            old = (self.p.coord, self.o5.coord,  tor1, tor2)
            # generate and check candidate positions
            j = 0
            while j < MAX_J:
                # construct P atom
                t1 = self.get_random_torsion(tor1,step)
                d = build_coord(self.c4v, self.c3v, self.o3v, O3_P_DIST, A1, t1, matrix1)
                self.p.coord = array((d[0], d[1], d[2]))
                self.pv = self.p.get_vector()
                matrix2 = get_ref_matrix(self.c3v, self.o3v, self.pv)
                k = 0
                while k < MAX_K:
                    # construct O5 atom
                    t2 = self.get_random_torsion(tor2,step)
                    d = build_coord(self.c3v, self.o3v, self.pv, P_O5_DIST, A2, t2, matrix2)
                    self.o5.coord = array((d[0], d[1], d[2]))
                    self.o5v = self.o5.get_vector()
                    # evaluate the candidate atoms
                    score = self.get_score(delta)
                    if score < delta:
                        delta = score
                        tor1, tor2 = t1, t2
                        j = MAX_J
                        k = MAX_K
                    k += 1
                j += 1
            if j == MAX_J:
                # old position is better
                self.p.coord, self.o5.coord, tor1, tor2 = old
                step = step*0.5
            i += 1
        

    def add_op12(self):
        """Adds OP1 OP2 atoms"""
        self.remove_op1_op2()
        phos = RNAChain('file', DATA_PATH+'phosphate_group.pdb')['1']
        #TODO: replace by Bio.PDB --> faster
        mobile = [phos['OP3'], phos['P'], phos["O5'"] ]
        fixed = [self.r1["O3'"], self.r2['P'], self.r2["O5'"] ]
        s = ModernaSuperimposer(fixed, mobile, phos.child_list)
        self.r2.add(phos['OP1'])
        self.r2.add(phos['OP2'])
    
    

class TerminalPhosphateBuilder(object):
    """
    Constructs missing P, OP1, OP2, OP3 atoms to a single residue.
    """
    def __init__(self,resi, structure):
        self.r1 = resi
        self.p = None
        self.o5 = self.r1["O5'"]
        self.o5v = self.o5.get_vector()
        self.c4 = self.r1["C4'"]
        self.c4v = self.c4.get_vector()
        self.c5 = self.r1["C5'"]
        self.c5v = self.c5.get_vector()
        self.c4 = self.r1["C4'"]
        self.c4v = self.c4.get_vector()
        self.ms = structure
    
    def build_phosphate_group(self):
        """Constructs a phosphate group based on O5', C5', and C4'"""
        self.remove_phosphate_group()
        atoms = self.__get_P_group_in_right_conformation()
        for at in atoms: 
            self.r1.add(at)

    def __has_local_clashes(self, ns, atoms):
        """Checks whether given atom list clashes with the neighbour atoms."""
        at_dict = {}
        for at in atoms:
            at_dict[at.name] = at

        if len(ns.search(at_dict['P'].get_coord(), 2.0, 'A')) > 1: return True # here we have connection with O5'
        if len(ns.search(at_dict['OP1'].get_coord(), 1.8, 'A')) > 0: return True
        if len(ns.search(at_dict['OP2'].get_coord(), 1.8, 'A')) > 0: return True
        if len(ns.search(at_dict['OP3'].get_coord(), 1.8, 'A')) > 0: return True
        return False

    def __get_P_group_in_right_conformation(self):
        """
        Checks possible B angles when building P group.
        Returns P group (4 atoms)
        that do not clash with rest of the structure.
        If there is no such conformation returns last tried conformation.
        """
        ns = NeighborSearch(self.ms.get_all_atoms()) 
        for ang in BETA:
            atoms = self.__get_phosphate_atoms(ang)
            if not self.__has_local_clashes(ns, atoms): 
                return atoms
        return atoms


    def __get_phosphate_atoms(self, beta):
        """
        Creates atoms (P, OP1, OP2, OP3) based on values from given resi and given B angle.
        Returns list of 4 created atoms.
        """
        p_coord = build_coord(self.c4v, self.c5v, self.o5v, P_O5_DIST, A3, beta)
        at_P = self.create_atom("P","P", p_coord.get_array())
        pv = at_P.get_vector()
        op1_coord = build_coord(self.c5v, self.o5v, pv, P_OP_DIST, P_OP_ANGLE, OP1_TORSION)
        op2_coord = build_coord(self.c5v, self.o5v, pv, P_OP_DIST, P_OP_ANGLE, OP2_TORSION)
        op3_coord = build_coord(self.c5v, self.o5v, pv, P_OP_DIST, P_OP_ANGLE, OP3_TORSION)
        at_OP1 = self.create_atom("OP1", "O", op1_coord.get_array())
        at_OP2 = self.create_atom("OP2", "O", op2_coord.get_array())
        at_OP3 = self.create_atom("OP3", "O", op3_coord.get_array())
        # might be usful here?
        # from Bio.PDB import calc_angle, calc_dihedral
        # angle = calc_angle(at1.get_vector(), at2.get_vector(), at3.get_vector())
        # dihedral = calc_dihedral(at1.get_vector(), at2.get_vector(), at3.get_vector(), at4.get_vector())
        return [at_P,  at_OP1,  at_OP2,  at_OP3]


    def create_atom(self, name, elem, coord):
        """Finds positions in between the O3' and C5'"""
        atom = Atom.Atom(name, coord, 0.0,1.0,' '," "+name,1, element=elem)
        return atom
        

    def remove_phosphate_group(self):
        """Removes eventually existing phosphate group atoms."""
        if self.r1.child_dict.has_key("P"):
            self.r1.detach_child("P")
        if self.r1.child_dict.has_key("OP1"):
            self.r1.detach_child("OP1")
        if self.r1.child_dict.has_key("OP2"):
            self.r1.detach_child("OP2")
        if self.r1.child_dict.has_key("OP3"):
            self.r1.detach_child("OP3")


