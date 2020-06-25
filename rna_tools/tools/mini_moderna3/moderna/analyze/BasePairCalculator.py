#!/usr/bin/env python
#
# BasePairCalculator.py
#
# Class for calculating interaction between RNA bases.
#
# http://iimcb.genesilico.pl/moderna/
#

"""
A module for  calculating the base pair type for two residues.
INPUT: two ModernaResidue objects
OUTPUT: a base pair type as a string, or None.
The main function here is base_pair_calc.
"""

__author__ = "Natalia Szostak"
__contributors__ = "Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Natalia Szostak"
__email__ = "ms_arwi@wp.pl"
__status__ = "alpha"


from rna_tools.tools.mini_moderna3.moderna.analyze.HBondCalculator import HBondCalculator
from rna_tools.tools.mini_moderna3.moderna.analyze.ChainConnectivity import is_ribose_complete
from numpy import array
from Bio.PDB.Atom import Vector

WC=["N1","N2","C2","N6","O6"]
S=["O2'","O3'","O4'","O5'","N9", "N3","C2","N2"]
H=["C8","N6","C5","C6","N7","O6"]
BB = ["OP1","OP2"]
EDGES_AG = {"WC":WC,"S":S,"H":H}


WC2=["C2","N3","N4","O2","O4"]
S2=["O2'","O3'","O4'","O5'","O2"]
H2=["N4","O4","C5","C6"]
EDGES_UC = {"WC":WC2,"S":S2,"H":H2}

NN_CUTOFF = 12.0

class BasePairInteraction(object):
    """Result from base pair calculation."""
    def __init__(self, resi1, resi2, bp_type):
        """Creates a base pair object."""
        self.resi1 = resi1
        self.resi2 = resi2
        self.type = bp_type
        
    def __repr__(self):
        return "%s %s %s"% \
            (self.resi1.identifier, self.type, self.resi2.identifier)
    
    @property
    def wobble(self):
        """True if a wobble pair."""
        if self.type == 'W/W':
            names = self.resi1.original_base, self.resi2.original_base
            if names in [('G', 'U'), ('U', 'G')]:
                return True
                
    @property
    def canonical(self):
        """True if a A-U or G-C pair (not wobble)."""
        return self.type in ['+/+', '-/-']

        
def get_edge_centers(resi):
    result = {}
    edges = resi.pyrimidine and EDGES_UC or EDGES_AG
    for edge in edges:
        n = 0
        coord = array([0.0, 0.0, 0.0])
        for name in edges[edge]:
            atom = resi.child_dict.get(name)
            if atom:
                coord += atom.coord
                n += 1
        if n>0:
            coord = coord/n
        result[edge] = Vector(coord)
    return result
        
def get_closest_edges(resi1, resi2):
    """Finds which edges are closest"""
    c1 = get_edge_centers(resi1)
    c2 = get_edge_centers(resi2)
    result = []
    for e1 in c1:
        for e2 in c2:
            distvec = c1[e1]-c2[e2]
            dist = distvec.norm()
            result.append((dist, e1, e2))
    result.sort()
    bp_type = result[0][1][0]+result[0][2][0]
    return bp_type

def validate_pair(resi1, resi2):
    """Checks if the two pairs are possibly paired."""
    if resi1.short_abbrev =='.' \
        or resi2.short_abbrev =='.'\
        or not is_ribose_complete(resi1) \
        or not is_ribose_complete(resi2):
            return False
    # skip residue pairs far away from each other
    nn_dist = resi1['N*'] - resi2['N*']
    if nn_dist > NN_CUTOFF:
        return False
    return True
    
def count_edges(atom, resi, ref_resi, edges, r1edge, r2edge):
    """Counts edges in h-bonds"""
    for edge in list(edges.keys()):
        if atom.get_name() in edges[edge]:
            if resi == ref_resi:
                r1edge.append(edge)
            else:
                r2edge.append(edge)

def __assign_edge(edge_list):
    """A function checking which edge dominates."""
    wc = edge_list.count("WC")
    s = edge_list.count("S")
    h = edge_list.count("H")
    
    # playing with > versus >= influences the  recall number a lot.
    if wc > s and wc > h:
        dominate = "W"
    elif s >= wc and s > h:
        dominate = "S"
    elif h > wc and h >= s:
        dominate = "H"
    else:
        dominate = "?"
    return dominate

def assign_bp_type(resi1, resi2, r1edge, r2edge, h_bond):
    """Determines base pair type from hbonds+edges."""
    # fix +/+ and -/- types
    hh = [str(hb) for hb in h_bond]
    
    resnames = resi1.original_base + resi2.original_base
    if resnames in ['GC', 'CG'] and 'N1-N3' in hh and 'N2-O2' in hh and 'N4-O6' in hh: 
        return '+/+' # G-C pairs
    if resnames in ['GU', 'UG'] and 'N1-O2' in hh and 'N3-O6' in hh: 
        return 'W/W' # G-U wobble pairs
    elif resnames in ['AU', 'UA'] and 'N3-N1' in hh and 'N6-O4' in hh: 
        return '-/-' # A-U pairs

    closest = get_closest_edges(resi1, resi2)    
    r1_type = __assign_edge(r1edge)
    if r1_type == '?': 
        r1_type = closest[0]
    r2_type = __assign_edge(r2edge)
    if r2_type == '?': 
        r2_type = closest[1]
    bp_type = r1_type + r2_type
    return bp_type

def base_pair_calc(resi1,resi2):
    """Assigns base pair type for two residues basing on HBonds calculated by HBondCalculator"""
    if not validate_pair(resi1, resi2):
        return None

    hb = HBondCalculator()
    h_bond = hb.calc_hbond_list(resi1, resi2) #calculating H-bonds
    if h_bond == []:# or len(h_bond) <1:
        return None
    else:
        r1edge = []
        r2edge = []
        for h in h_bond:
            edges = h.donor_resi.pyrimidine and EDGES_UC or EDGES_AG
            count_edges(h.donor, h.donor_resi, resi1, edges, r1edge, r2edge)
            edges = h.acceptor_resi.pyrimidine and EDGES_UC or EDGES_AG
            count_edges(h.acceptor, h.acceptor_resi, resi1, edges, r1edge, r2edge)
            
        bp_type = assign_bp_type(resi1, resi2, r1edge, r2edge, h_bond)
        return BasePairInteraction(resi1, resi2, bp_type)


