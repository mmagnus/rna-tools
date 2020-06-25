#!/usr/bin/env python
#
# mol_parser.py
#
# Reads molecule structures into a data structure
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, Kristian Rother"
__credits__ = ["Sabrina Hofmann"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

from .MolParameters import *
import re
    

class Bond:
    """Something connecting two atoms"""
    def __init__(self,atom1,atom2,valence, special=(0,0,0)):
        self.atom1 = atom1
        self.atom2 = atom2
        self.valence = valence

class Atom(dict):
    """Atom that has several default and lots of optional properties."""
    bond_class = Bond
    
    def __init__(self, element):
        """
        element - string indicating an element from the periodic system.
        """
        dict.__init__(self)
        self.element = element
        self.bonds = []
        self.n_bonds = 0
        self.added_valences = 0
        
        self['atom_name']  = element
        self['attributes'] = []
        
        self._is_fully_bonded = False # SPEEDUP
    
    def add_bond(self, atom, valence, special=(0,0,0)):
        self.bonds.append(self.bond_class(self,atom, valence, special))
        self.added_valences += valence
        self.n_bonds += 1
        
    def add_attribute(self, att):
        if att not in self['attributes']:
            self['attributes'].append(att)
            
    def is_fully_bonded(self):
        if self._is_fully_bonded: return True
        elif self.element=='O' and (len(self.bonds)>=2 or self.added_valences>=4): self._is_fully_bonded = True # because of aromatic bonds
        elif self.element=='C' and (len(self.bonds)>=4 or self.added_valences>=6): self._is_fully_bonded = True
        elif self.element=='N' and (len(self.bonds)>=4 or self.added_valences>=6): self._is_fully_bonded = True
        return self._is_fully_bonded
        
    def __str__(self):
        return self.get('bondschema',self.element)
        
    def detect_bonding_schema(self):
        """
        creates a string that summarizes number and element of neighbors
        e.g. 'CC1H111'
        and the ProtOr atomtype according to Tsai 1999
        """
        neighbors = {}
        for bond in self.bonds:
            elem = bond.atom2.element
            neighbors.setdefault(elem,[])
            neighbors[elem].append(str(bond.valence))
        k = list(neighbors.keys())
        k.sort()
        atomtype = self.element
        for kk in k:
            # make a string of the form CC1O12
            neighbors[kk].sort()
            atomtype += kk + ''.join(neighbors[kk])
        self['bondschema'] = atomtype

        protortype = self.element
        if atomtype[0] != 'H':
            bonded_h = 0
            if 'H' in neighbors: bonded_h = len(neighbors['H'])
            protortype += "%iH%i"%(len(self.bonds),bonded_h)

    def get_molstring(self, taboo_list,depth=-1):
        """
        Recursive procedure that builds a string from bond schema strings.
        The entire sub-molecule appears there,
        e.g. C1111(C1111(H1,H1,H1),H1,H1) for an ethyl group.
        Requires the 'bondschema' key to be set for full information.
        The parts of the string are sorted according to alphabet.
        b) length.
        """
        if depth == 0: return ""
        new_taboo = taboo_list[:]
        
        s = str(self)
        neighbors = []
        for bond in self.bonds:
            # update list of visited atoms
            if bond.atom2 not in taboo_list: 
                new_taboo.append(bond.atom2)
        for bond in self.bonds:
            if bond.atom2 not in taboo_list: 
                neighbors.append(bond.atom2.get_molstring(new_taboo,depth-1))
        neighbors.sort()
            
        if neighbors:
            s += "("+','.join(neighbors)+')'
        return s


class Molecule(list):
    """
    Contains a data structure with atoms and bonds that can be created from
    a .mol file, a .mol2 file or from a Bio.PDB.Residue object.
    """
    atom_class = Atom

    def __init__(self):
        list.__init__(self)
        self.conjugated_systems = []

    def convert_pdb_atom(self, atom):
        """Returns a Atom object from a PDB.Atom object"""
        element = atom.id[0]
        if element in "1234567890": element = atom.id[1]
        if len(atom.id)>1 and element == 'S' and atom.id[1].upper() == 'E':element = 'Se'
        at = self.atom_class(element)
        if element == 'H': 
            return None
        at['atom_name'] = atom.id
        return at

    def parse_resi(self,resi):
        """
        Creates a molecule object from a Bio.PDB.Residue object.
        Crude solution that has to be checked for better distance constraints.
        Do not use it for other things as nucleotide bases.
        Hydrogens are skipped.
        """
        atoms = [atom for atom in resi if not atom.name.strip()[0] =='H']
        # create atoms
        for atom in atoms:
            at = self.convert_pdb_atom(atom)
            if at: self.append(at)

        # create bonds
        for i, atom1 in enumerate(self):
            if atom1.is_fully_bonded(): continue
            for j in range(i+1, len(self)):
                atom2 = self[j]
                if atom2.is_fully_bonded(): continue
                distance= atoms[i]-atoms[j]
                if distance < MAX_BOND_DISTANCE:
                    min_sgl, max_sgl = SINGLE_BOND_DISTANCES.get(atom1.element,{}).get(atom2.element,DEFAULT_BOND_LENGTHS[0])
                    min_dbl, max_dbl = DOUBLE_BOND_DISTANCES.get(atom1.element,{}).get(atom2.element,DEFAULT_BOND_LENGTHS[1])
                    min_tri, max_tri = TRIPLE_BOND_DISTANCES.get(atom1.element,{}).get(atom2.element,DEFAULT_BOND_LENGTHS[2])
                    valence = 0
                    if min_sgl < distance < max_sgl:
                        valence = 1 # single bond
                    elif min_dbl < distance < max_dbl:
                        valence = 2 # double bond
                    elif min_tri < distance < max_tri:
                        valence = 3 # triple bond
                    if valence:
                        atom1.add_bond(atom2,valence,(0,0,distance))
                        atom2.add_bond(atom1,valence,(0,0,distance))
        
        # print debug information
        if 0:
            for atom1, atom2 in zip(atoms, self):
                atom2.detect_bonding_schema()
                print((atom1.fullname, atom2['bondschema']))

    def parse_molfile(self,filename):
        """
        Reads a .mol file and parse the contents.
        """
        for l in open(filename):
            if re.search("^\s+\-*\d+\.\d+",l):
                # parse line with atom coordinates
                t = re.sub("\s+","\t",l[:-1]).split('\t')
                atom = self.atom_class(element=t[4])
                self.append(atom)
                atom['coordinates'] = (float(t[1]),float(t[2]),float(t[3]))
                atom['rdf_index'] = int(l[60:63]) # atom mapping from BioPath RDF files
                
            elif re.search("^\s*\d+\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+(\s+\d+)*\s*\n*\Z",l):
                # parse bonds line
                atom1 = self[int(l[0:3])-1]
                atom2 = self[int(l[3:6])-1]
                valence = int(l[6:9])
                special = (int(l[9:12]),int(l[12:15]),int(l[15:18]))
                # the special features indicate the orientation of the sugar hydroxyl groups.
                atom1.add_bond(atom2,valence,special)
                atom2.add_bond(atom1,valence,special)
            elif re.search('^...END',l): 
                return
                
    def parse_mol2(self,filename):
        """
        Read a .mol2 file and parse the contents.
        """
        for l in open(filename):
            if re.search("^[\s\d]{7}\s[^\s]+\d+",l):
                # parse coordinates line        
                t = string.split(re.sub("\s+","\t",l[:-1]),'\t')
                atom = self.atom_class(element=re.sub('\d','',t[2]))
                self.append(atom)
                atom['coordinates'] = (float(l[17:26]),float(l[27:36]),float(l[37:46]))
                atom['rdf_index'] = int(int(l[:7])) # atom mapping from BioPath RDF files
                
            elif re.search("^\s*\d+\s+\d+\s+\d+\s+[^\s]+\s*\n*\Z",l):
                # parse bonds line
                atom1 = self[int(l[6:12])-1]
                atom2 = self[int(l[12:18])-1]
                valence = l[21:23]
                if valence == 'ar':
                    valence=2
                    atom1['attributes'].append('aromatic')
                    atom2['attributes'].append('aromatic')
                elif valence == 'am':
                    valence=2
                else:
                    valence = int(valence)
                atom1.add_bond(atom2,valence)
                atom2.add_bond(atom1,valencel)
                


