#!/usr/bin/env python
#
# GeometryStatistics.py
#
# Statitic for base type recognition.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Kristian Rother"
__copyright__ = "Copyright 2008, Genesilico"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "beta"

from Bio.PDB.vectors import Vector, calc_angle, calc_dihedral
import os, math, re
    
class GeometryError(Exception): pass

class AtomDefinition(object):
    """
    Contains a single atom definition, like
    O2'
    X:O2'
    X[G]:N1
    X[A]:N1
    X{A}:N1
    X+1{A}[G]:N1
    NUC:Y+1[C]:O2'
    Returns an AtomDefinition object.
    """
    def __init__(self,atom_def):
        self.atom_name = None
        self.res_code = None
        self.res_name = None
        self.res_offset = 0
        self.chain_id = None
        self.mol_type = None
        self.parse_atom_definition(atom_def)
    
    def parse_atom_definition(self,atom_def):
        # first element may be omitted
        token = atom_def.split(':')
        if not 2 <= len(token) <= 3:
            raise GeometryError("invalid atom definition: '%s'"%atom_def)
            
        self.atom_name = token[-1] # atom name is always the last
        self.res_code = token[-2][0] # residue comes second last
        suffix = token[-2][1:]
        if len(suffix)>1:
            #chain = re.findall('\{(.+)\}',suffix)
            resi = re.findall('\[(.+)\]',suffix)
            offset = re.findall('^([+-]\d+)',suffix)
            #if chain: self.chain_id = chain[0]
            if resi: self.res_name = resi[0]
            if offset: self.res_offset = int(offset[0]) 
        if len(token) == 3:
            self.mol_type = token[0]
            
    def get_atom(self,residue_list, index):
        """Returns a single atom"""
        ofs_index = index + self.res_offset
        if 0 <= ofs_index < len(residue_list):
            # assert ofs_index is still in the same change and resno matches.
            pass
            # get the residue
            resi = residue_list[ofs_index]
            # assert the chain is set correctly
            if not self.chain_id:# or resi.parent.id == self.chain_id:
                # assert the residue name is set correctly
                resname_ok = True
                if self.res_name:
                    rn = resi.resname.strip()
                    if self.res_name!=rn: resname_ok = False
                    if self.res_name == 'R' and rn in ['A','G'] \
                       or self.res_name == 'Y' and rn in ['C','U']:
                        resname_ok = True
                if resname_ok:
                    # get the atom
                    if resi.has_id(self.atom_name):
                        return resi[self.atom_name]
        return None

            
class GeometryExpression(list):
    """Handles the "X:C5',X:C4'"-like strings given to the GeometryStatistics class. 
    Is a list of AtomDefinition objects."""
    
    def __init__(self,expression):
        """Parses the expression and creates AtomDefinition objects."""
        list.__init__(self)
        for atom_def in expression.split(','):
            self.append(AtomDefinition(atom_def))
            
        self.residue_codes = self.get_residue_codes()
        
    def get_residue_codes(self):
        """Collects all residue codes occurring in this expression."""
        codes = []
        for ad in self:
            if ad.res_code not in codes:
                codes.append(ad.res_code)
        return codes
            
    def get_residues_by_codes(self, struc):
        """Returns a nested list of residues, checking the chain type of each."""
        residues = []
        for rc in self.residue_codes:
            resi_list = [resi for resi in struc]
            residues.append(resi_list)
        return residues
        
    def has_distinct_residues(self, residues, r_index):
        """Returns true only if the indexed residues are not identical."""
        # MR added this. KR please check
        # checks whwther there are any residues
        if len(residues)==0 or (len(residues)==1 and residues[0]==[]): return False
        
        # get the indexed residues
        resi = []
        for i in range(len(r_index)):
            resi.append(residues[i][r_index[i]])
        # look for duplicates
        for i in range(len(resi)):
            if resi[i] in resi[i+1:]:
                    return False
        return True
    
    def get_atom_combination(self, residues, r_index):
        """Returns a list of atoms corresponding to the list of atom definitions.
        residues is a nested list of residue objects,
        r_index indicates the positions of each of the nested residues.
        """
        atoms = []
        for ad in self:
            # find the right residue and get an atom for it.
            column = self.residue_codes.index(ad.res_code)
            resi_list = residues[column]
            atoms.append(ad.get_atom(resi_list,r_index[column]))
        return atoms
        
    def check_atom_combination(self,atoms):
        if None in atoms: 
            return False
        return True

    def get_atoms(self, struc):
        residues = self.get_residues_by_codes(struc)
        r_index = [0] * len(residues)
        ready = False
        while not ready:
            if self.has_distinct_residues(residues,r_index): 
                # retrieve atom information and generate if it is ok.
                atoms = self.get_atom_combination(residues, r_index)
                if self.check_atom_combination(atoms):
                    yield atoms
            # count up atom indices
            i = 0
            while i< len(r_index):
                r_index[i] += 1
                if r_index[i] >= len(residues[i]):
                    r_index[i] = 0
                    i += 1
                    if i == len(r_index):
                        ready = True
                else:
                    i = len(r_index)
                    

class GeometryResult(list):
    """List of floats producing tabular reports and plots."""
    def __init__(self, title, angles=False):
        """
        title - string
        angles - if set true, 360.0 = 0.0, taken into account in all average/stddev calculations.
        """
        list.__init__(self)
        self.title = title
        self.angles = angles
        
    def __str__(self):
        result = "%s\t%i\t%7.3f\t%7.3f\n"%(self.title, len(self), self.get_avg(), self.get_stddev())
        return result
        
    def append(self, value):
        if self.angles :
            assert 0.0 <= value[0] < 360.0
        list.append(self,value)
        
    def get_avg(self):
        """returns the average value"""
        avg = sum([d[0] for d in self])/len(self)
        if self.angles:
            # angle overlap correction
            dev = sum([abs(d[0]-avg) for d in self])
            alt_avg = avg<180.0 and avg + 180.0 or avg-180.0
            alt_dev = sum([min(abs(d[0]-alt_avg),abs(d[0]-360.0-alt_avg)) for d in self])
            avg = dev<alt_dev and avg or alt_avg
        return avg
        
    def get_stddev(self):
        """returns the standard deviation."""
        avg = self.get_avg()
        if len(self) > 1:
            if self.angles:
                sq_sum = 0.0
                for d in self:
                    diff = min(abs(d[0]-avg), abs(d[0]-360.0-avg))**2
                    sq_sum += diff
            else:
                sq_sum = sum([(d[0]-avg)**2 for d in self])
            stddev = math.sqrt(sq_sum/(len(self)-1.0))
        else:
            stddev = -1.0
        return stddev
        
    def add_data(self,data):
        for d in data:
            self.append(d)
       
    def write_plot(self, fn, xmin=0.0, ymin=0.0, xmax=None, ymax=None):
        pass
        
    def write_table(self, fn):
        out = [self.title + '\n']
        for d in self:
            line = "%8.3f\t%s\n"%(d[0],'\t'.join(d[1:]))
            out.append(line)
        open(fn,'w').writelines(out)
        
    def histogram(self):
        pass
        
    def scatterplot(self):
        pass


class GeometryStatistics(object):
    """
    Calculates geometry statistics for a single structure object.
    """
    def __init__(self, structure):
        """structure is a ModernaStructure object."""
        self.struct = structure
        
    def get_structures(self):
        yield self.struct
                
    def get_atom_annotation(self,atom):
        res = atom.parent
        return '%s;'%(res.identifier)
        #rid = str(res.id[1])+res.id[2].strip()
        #chain = res.parent
        #cid = chain.id
        #struc = chain.parent.parent.id
        #anno = '%s/%s/%s;'%(struc,cid,rid)
        #return anno
                
    def get_distances(self, expression):
        result = GeometryResult('distances of (%s)'%expression)
        expr = GeometryExpression(expression)
        for struc in self.get_structures():
            for atom1,atom2 in expr.get_atoms(struc):
                a1 = self.get_atom_annotation(atom1)
                a2 = self.get_atom_annotation(atom2)
                result.append((atom1-atom2,a1,a2))
        return result
        
    def get_angles(self, expression):
        result = GeometryResult('angles of (%s)'%expression, angles=True)
        expr = GeometryExpression(expression)
        for struc in self.get_structures():
            for atom1,atom2,atom3 in expr.get_atoms(struc):
                a = calc_angle(atom1.get_vector(),atom2.get_vector(),atom3.get_vector())
                #a=angle(atom2.coord-atom1.coord,atom2.coord-atom3.coord)
                a = a*180.0/math.pi
                a1 = self.get_atom_annotation(atom1)
                a2 = self.get_atom_annotation(atom2)
                a3 = self.get_atom_annotation(atom3)
                result.append((a,a1,a2,a3))
        return result
    
    def calc_dihedral(self,atom1,atom2,atom3,atom4):
        dihedral = calc_dihedral(atom1.get_vector(),\
            atom2.get_vector(),atom3.get_vector(), atom4.get_vector())
        dihedral *= 180/math.pi
        if dihedral < 0: dihedral += 360
        return dihedral
        
    def get_dihedrals(self, expression):
        result = GeometryResult('dihedrals of (%s)'%expression, angles=True)
        expr = GeometryExpression(expression)
        for struc in self.get_structures():
            for atom1,atom2,atom3,atom4 in expr.get_atoms(struc):
                try:
                    a1 = self.get_atom_annotation(atom1)
                    a2 = self.get_atom_annotation(atom2)
                    a3 = self.get_atom_annotation(atom3)
                    a4 = self.get_atom_annotation(atom4)
                    result.append((self.calc_dihedral(atom1,atom2,atom3,atom4),a1,a2,a3,a4))
                except ValueError as e:
                    print((str(e)))
                    print((atom1.coord,atom2.coord,atom3.coord,atom4.coord))
        return result
        
class PDBSetGeometryStatistics(GeometryStatistics):
    """
    Calculates statistics over a set of structure files.
    """
    def __init__(self, pdb_path):
        if pdb_path[-1] != os.sep: 
            pdb_path += os.sep # make sure there is a terminating slash
        self.pdb_path = pdb_path
        
    def get_structures(self):
        from moderna import ModernaStructure
        for fn in os.listdir(self.pdb_path):
            if fn[-4:].upper() in ['.ENT','.PDB']:
                if fn[-6] == '_': chain_id = fn[-5]
                else: chain_id = 'A'
                structure=ModernaStructure('file', self.pdb_path + fn, chain_name=chain_id)
                yield structure
    

if __name__ == '__main__':
    print("""
    
GeometryStatistics by Kristian.
(c) 2008 Genesilico

alpha version (!)

Calculates statistics for atoms in PDB files.

works like this:

from geometry_statistics import GeometryStatistics
gs = GeometryStatistics('path_with_pdb_files')
result = gs.get_distances("X:C4',X:C1'")
result.write_table('out.txt')
result.write_plot('out.png')

also implemented:
result = gs.get_angles("X:C4',X:C3',X:C1'")
result = gs.get_dihedrals("X:C4',X:C3',X:O3',X:C1'")

there are a number of ways how atoms can be specified. Examples:
   O2'
    X:O2'
    X[G]:N1
    X[A]:N1
    Y+1:O2'
    NUC:Y+1[C]:O2'

Things to be checked:
- Are the angles calculated properly (between vector(atom1-atom2) and vector(atom1-atom3)?

Known issues:
- calculates junk for multi-chain structures
- cannot distinguish between protein/nucleic acid chains -> slower than it could be
- plotting does not work yet.

Annotation of atoms:
In the result tables, each atom used for calculation in each line is annotated
by a string separable by ';' like this:

1EHZ/A/25;C3'-endo;<<(A/24);cis/+/+(A/10);

The first segment '1EHZ/A/25', describes the residue where the atom is from.
The second segment 'C3'-endo', defines the pucker of the ribose.
The third segment '<<(A/24)', describes a stacking interaction.
The fourth segment 'cis/+/+(A/10)' describes a base pair.

If no annotation was found, some of these descriptors may be missing.
Also, a single residue may undergo more than one stacking/base pair interaction.

""")
