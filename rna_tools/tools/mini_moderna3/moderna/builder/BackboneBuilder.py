#!/usr/bin/env python
#
# BackboneBuilder.py
#
# For backbone fixing.
#
# http://iimcb.genesilico.pl/moderna/
#

"""
BackboneBuilder.py

module for reconstructing the full backbone 
between two riboses.
e.g. after loop insertsions.
"""

__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"


from rna_tools.tools.mini_moderna3.moderna.analyze.GeometryParameters import GeometryStandards
from rna_tools.tools.mini_moderna3.moderna.builder.PhosphateBuilder import PhosphateBuilder
from rna_tools.tools.mini_moderna3.moderna.RNAResidue import RNAResidue
from rna_tools.tools.mini_moderna3.moderna.builder.FCCDLoopCloser import FCCDLoopCloser
from rna_tools.tools.mini_moderna3.moderna.builder.CoordBuilder import build_coord
from rna_tools.tools.mini_moderna3.moderna.analyze.RNASuites import TORSIONS, DEFAULT_TORSIONS
from rna_tools.tools.mini_moderna3.moderna.analyze.ChainConnectivity import are_residues_connected, \
    is_backbone_intact, is_phosphate_intact, is_backbone_congested
from numpy import array
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log

GEO_STD = GeometryStandards()
O3_P_DIST  = GEO_STD.get_standard("X:O3',X+1:P")
P_O5_DIST  = GEO_STD.get_standard("X:P,X:O5'")
O5_C5_DIST = GEO_STD.get_standard("X:O5',X:C5'")
C5_C4_DIST = GEO_STD.get_standard("X:C5',X:C4'")
C4_C3_DIST = GEO_STD.get_standard("X:C4',X:C3'")
C3_C2_DIST = GEO_STD.get_standard("X:C3',X:C2'")


class BackboneBuilder(object):
    """
    combines various strategies to rebuild a backbone.
    """
    def __init__(self, resi1, resi2, struc):
        self.resi1 = resi1
        self.resi2 = resi2
        self.struc = struc
        if not self.is_intact():
            log.write_message('''The backbone between residues %s and %s is not\
 intact. Attempting rebuilding.'''%(resi1.identifier, resi2.identifier))
            self.build()
            
    def write_bb_status(self):
        log.write_message(str(self.is_intact()))
        self.write_resi_status(self.resi1)
        self.write_resi_status(self.resi2)
        if are_residues_connected(self.resi1, self.resi2):
            log.write_message('connection between both residues OK')
        else:
            log.write_message('connection between both residues BROKEN')
        
    def write_resi_status(self, resi):
        log.write_message('Status of residue: %s'%resi.identifier)
        if is_backbone_intact(resi):
            log.write_message('\tbackbone  OK')
        else:
            log.write_message('\tbackbone  INTERRUPTED')
        if is_phosphate_intact(resi):
            log.write_message('\tphosphate OK')
        else: log.write_message('\tphosphate BROKEN')
        if is_backbone_congested(resi):
            log.write_message('\tclashes   OCCUR')
        else:
            log.write_message('\clashes   OK')
    
    def is_intact(self):
        """Checks whether the backbone is complete."""
        if not is_backbone_intact(self.resi1, mode="3'") or \
            not is_backbone_intact(self.resi2, mode="5'") or \
            not are_residues_connected(self.resi1, self.resi2) or \
            not is_phosphate_intact(self.resi2) or \
            is_backbone_congested(self.resi1, mode="3'") or \
            is_backbone_congested(self.resi2, mode="5'"):
                return False
        return True

    def create_resi_from_torsions(self, torsions):
        """Creates residue with backbone atoms in different positions."""
        resi2_copy = RNAResidue(self.resi2, \
                alphabet_entry=self.resi2.alphabet_entry)
        epsilon, zeta, alpha, beta, gamma, delta = torsions
        c4c3 = self.resi2["C4'"]-self.resi2["C3'"]
        c3c2 = self.resi2["C3'"]-self.resi2["C2'"]
        at_p = build_coord(self.resi1["C4'"].get_vector(), \
                self.resi1["C3'"].get_vector(), \
                self.resi1["O3'"].get_vector(), O3_P_DIST, 109.0,epsilon)
        at_o5 = build_coord(self.resi1["C3'"].get_vector(), \
                self.resi1["O3'"].get_vector(), \
                at_p, P_O5_DIST, 109.0, zeta)
        at_c5 = build_coord(self.resi1["O3'"].get_vector(), \
                at_p, at_o5, O5_C5_DIST, 109.0, alpha)
        at_c4 = build_coord(at_p, at_o5, at_c5, C5_C4_DIST, 109.0, beta)
        at_c3 = build_coord(at_o5, at_c5, at_c4, c4c3, 109.0, gamma)
        at_c2 = build_coord(at_c5, at_c4, at_c3, c3c2, 109.0, delta)
        resi2_copy['P'].coord = array((at_p[0], at_p[1], at_p[2]))
        resi2_copy["O5'"].coord = array((at_o5[0], at_o5[1], at_o5[2]))
        resi2_copy["C5'"].coord = array((at_c5[0], at_c5[1], at_c5[2]))
        resi2_copy["C4'"].coord = array((at_c4[0], at_c4[1], at_c4[2]))
        resi2_copy["C3'"].coord = array((at_c3[0], at_c3[1], at_c3[2]))
        resi2_copy["C2'"].coord = array((at_c2[0], at_c2[1], at_c2[2]))
        return resi2_copy
        
    def run_fccd(self, torsions=DEFAULT_TORSIONS, threshold=0.05, maxit=200):
        """Runs the FCCDLoopCloser algorithm."""
        r1 = self.resi1
        r2 = self.resi2
        r2a = self.create_resi_from_torsions(torsions)
        fixed = [r2["C4'"], r2["C3'"], r2["C2'"]]
        moving = [
                  r1["C4'"], r1["C3'"], r1["O3'"], 
                  r2a["P"], r2a["O5'"], r2a["C5'"], r2a["C4'"], 
                  r2a["C3'"], r2a["C2'"]
                  ]
        lc = FCCDLoopCloser(moving, fixed)
        msg, rmsd, it = lc.run_fccd(threshold, maxit)
        lc.copy_vectors_to_atoms()
        # edit atoms that dont move
        r2a.detach_child("C4'")
        r2a.detach_child("C3'")
        r2a.detach_child("C2'")
        r2a.add(r2["C4'"])
        r2a.add(r2["C3'"])
        r2a.add(r2["C2'"])
        # add moved residue to structure
        self.struc.remove_residue(r2.identifier)
        self.struc.add_residue(r2a)
        self.resi2 = self.struc[r2a.identifier] # gets clone of r2a!
        log.write_message('.. result: %s'%msg)
        log.write_message('     rmsd: %6.3f   iterations: %i'%(rmsd, it))
        
    def build_phosphate(self):
        """tries to remodel the phosphate."""
        log.write_message('.. remodeling phosphate.')
        phosb = PhosphateBuilder(self.resi1, self.resi2)
        phosb.build()
        
    def build_fccd(self):
        """Tries to remodel a suite without O3'."""
        log.write_message(".. remodeling suite between O3' and C4'.")
        self.run_fccd()
        phosb = PhosphateBuilder(self.resi1, self.resi2)
        phosb.add_op12()

    def build_suites(self):
        """Tries different combinations of starting angles."""
        log.write_message(".. trying %i different suite torsion combinations."\
            %(len(TORSIONS)-1))
        for torsions in TORSIONS[1:]:
            self.run_fccd(torsions=torsions, threshold=0.2, \
                maxit=100)
            phosb = PhosphateBuilder(self.resi1, self.resi2)
            phosb.add_op12()
            if self.is_intact(): return

    def build(self):
        """Runs procedures to construct backbone between two residues."""
        methods = [self.build_phosphate, self.build_fccd, self.build_suites]
        for method in methods:
            method()
            if self.is_intact():
                log.write_message('Backbone reconstruction successful.')
                return
        
        self.write_bb_status()
        log.write_message('WARNING: Backbone reconstruction failed. \
The structure needs to be refined.\n')

