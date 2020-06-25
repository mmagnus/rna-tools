#!/usr/bin/env python
#
# GeometryAnanlyzer.py
#
# Analyze geometry of a RNA structure.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, Genesilico"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"


from numpy import array
from rna_tools.tools.mini_moderna3.moderna.analyze.GeometryStatistics import AtomDefinition, GeometryStatistics
from rna_tools.tools.mini_moderna3.moderna.analyze.GeometryParameters import GeometryStandards
              
class GeometryAnalyzer:
    standards = GeometryStandards()
    
    def __init__(self, structure):
        """Takes a ModernaStructure object"""
        self.struc = structure
        self.bad_bonds = []
        self.bad_angles = []
        self.bad_dihedrals = []
    
    def analyze(self):
        gs = GeometryStatistics(self.struc)
        self.check_bonds()
        self.check_angles()
        self.check_dihedrals()
        
                
    def find_outliers(self, descriptor, values):
        result = []
        for v in values:
            if self.standards.is_outlier(descriptor, v[0]):
                result.append((descriptor, v))
        return result

    def check_bonds(self):
        gs = GeometryStatistics(self.struc)
        self.bad_bonds = []
        for b in self.standards.bonds:
            result = gs.get_distances(b)
            self.bad_bonds += self.find_outliers(b, result)
        return self.bad_bonds
        
    def check_angles(self):
        gs = GeometryStatistics(self.struc)
        self.bad_angles = []
        for a in self.standards.angles:
            result = gs.get_angles(a)            
            self.bad_angles += self.find_outliers(a, result)
        return self.bad_angles
        
    def check_dihedrals(self):
        gs = GeometryStatistics(self.struc)
        self.bad_dihedrals = []
        for d in self.standards.dihedrals:
            result = gs.get_dihedrals(d)
            self.bad_dihedrals += self.find_outliers(d, result)
        return self.bad_dihedrals
        
    def get_resi_num_str(self,b):
        """ """
        resi = set(b[1][1:])
        if len(resi)==1:
            return '- Residue %s'%((b[1][1]).replace(';',''))
        else:
            res_num = list(resi)
            return '- Residues %s --- %s'%(res_num[0].replace(';',''),res_num[1].replace(';',''))

    def __str__(self):
        result = ""
        if self.bad_bonds: result += 'Unusual bond lengths:\n'
        for b in self.bad_bonds:
            result += '%s (%s): %.02f\n' %(self.get_resi_num_str(b), b[0], b[1][0])
        if self.bad_angles: result += "\nUnusual bond angles:\n"
        for b in self.bad_angles:
            result += '%s (%s): %.02f\n' %(self.get_resi_num_str(b), b[0], b[1][0])
        if self.bad_dihedrals: result += "\nUnusual dihedral angles:\n"
        for b in self.bad_dihedrals:
            result +=  '%s (%s): %.02f\n' %(self.get_resi_num_str(b), b[0], b[1][0])
        return result+'\n'
        
    def __repr__(self):
        return str(self)
        
        
