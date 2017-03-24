from numpy import array, sqrt

from Bio import PDB
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO, Superimposer

#from Bio.SVDSuperimposer import SVDSuperimposer
#from SVDSuperimposerDontMove import SVDSuperimposer
import numpy

class BasePair:
    """BasePair Class

    methods:
    - get_coord()
    - __sub__() - calc distance between two base pairs

    attributs:
    - coord - centroid of base pair"""
    def __init__(self, a,b,struct=None):
        """a,b, struct biopython residue

        struct is needed for saving a pdb file (default is None)"""
        self.a = a
        self.b = b
        self.struct = struct
        self.coord = None
        self.name = 'bp-' + a.get_resname().strip() + str(a.get_id()[1]) + '-' + b.get_resname().strip() + str(b.get_id()[1])
        self.coord = self.calc_coord()
        self.atoms_ab = self.get_atoms_ab()
        self.atoms_ba = self.get_atoms_ba()

    def get_name(self):
        return self.name
    def calc_coord(self):
        """Return coord of base pair

        exclude sugar and posphate and hydrogens"""
        centroid=PDB.Vector(0,0,0)
        counter=0.0
        #exclude = ["O5'","C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "P", "OP2", "OP1", "OP3","C1'"] # 
        exclude = []
        for atom in self.a:
            if atom.name not in exclude:
                if not atom.name.startswith('H'): 
                    centroid=centroid+atom.get_vector()
                    counter+=1
        for atom in self.b:
            if atom.name not in exclude:
               if not atom.name.startswith('H'): 
                    centroid=centroid+atom.get_vector()
                    counter+=1

        centroid=centroid/counter
        self.coord = numpy.array((tuple(centroid)), 'f')
        return self.coord

    def __sub__(self, other):
        """Calc distance between two base pairs"""
        diff = self.coord - other.coord 
        return numpy.sqrt(numpy.dot(diff, diff)) 
    def __str__(self):
        return self.name

    def overlap_with(self, other):
        diff = self.coord - other.coord         
        self.new_coord = self.coord - diff

        self.a_bk = self.a.copy()
        self.b_bk = self.b.copy()

        for atom in self.a:
            atom.set_coord(atom.get_coord() - diff)
        for atom in self.b:
            atom.set_coord(atom.get_coord() - diff)

    def save(self):
        if not self.struct:
            raise Exception('self.struct was not defined! Can not save a pdb!')
        class BpSelect(Select):
            def accept_residue(self, residue):
                if residue.get_id()[1] == 1 or residue.get_id()[1] == 43:
                    return 1
                else:
                    return 0

        io = PDBIO()
        io.set_structure(self.struct)
        fn = self.name + '.pdb'
        io.save(fn, BpSelect())
        return 'Saved to: %s ' % fn

    def get_atoms_ab(self):
        """Get 18 atoms per base pair.

        18
        [<Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C1'>, <Atom C2'>, <Atom O2'>, <Atom N9>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C1'>, <Atom C2'>, <Atom O2'>, <Atom N1>]"""
        self.atoms = []

        #for atom in self.a:
        #    atoms.append(atom)
        #for atom in self.a:
        #    atoms.append(atom) 
        include = ["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]# "OP2", "OP1", "OP3",  "P" "O5'" # <- first pair does not have them
        #self.atoms.extend([self.a['C3\''], self.b['C3\'']])

        for at in self.b:
            if at.name in include:
                self.atoms.append(at)
            if self.b.resname.strip() in ['G', 'A'] and at.name == "N9":
                self.atoms.append(at)
            if self.b.resname.strip() in ['U', 'C'] and at.name == "N1":
                self.atoms.append(at)

        for at in self.a:
            if at.name in include:
                self.atoms.append(at)
            if self.a.resname.strip() in ['G', 'A'] and at.name == "N9":
                self.atoms.append(at)
            if self.a.resname.strip() in ['U', 'C'] and at.name == "N1":
                self.atoms.append(at)

        return self.atoms

    def get_atoms_ba(self):
        """Return 18 atoms per base pair.

        18
        [<Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C1'>, <Atom C2'>, <Atom O2'>, <Atom N9>, <Atom C5'>, <Atom C4'>, <Atom O4'>, <Atom C3'>, <Atom O3'>, <Atom C1'>, <Atom C2'>, <Atom O2'>, <Atom N1>]

        print len(coords1) --> 18!"""
        self.atoms = []

        #for atom in self.a:
        #    atoms.append(atom)
        #for atom in self.a:
        #    atoms.append(atom) 
        include = ["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]# "OP2", "OP1", "OP3",  "P" "O5'" # <- first pair does not have them
        #self.atoms.extend([self.a['C3\''], self.b['C3\'']])

        for at in self.a:
            if at.name in include:
                self.atoms.append(at)
            if self.a.resname.strip() in ['G', 'A'] and at.name == "N9":
                self.atoms.append(at)
            if self.a.resname.strip() in ['U', 'C'] and at.name == "N1":
                self.atoms.append(at)

        for at in self.b:
            if at.name in include:
                self.atoms.append(at)
            if self.b.resname.strip() in ['G', 'A'] and at.name == "N9":
                self.atoms.append(at)
            if self.b.resname.strip() in ['U', 'C'] and at.name == "N1":
                self.atoms.append(at)

        return self.atoms

    def calc_rmsd_to(self, other):
        """Return rmsd to other bp (it calc two rmsd a-b and b-a and return min)."""
        #r1
        coords1 = array(tuple([tuple(atom.get_vector()) for atom in self.atoms_ab]), 'f')
        coords2 = array(tuple([tuple(atom.get_vector()) for atom in other.atoms_ab]), 'f')

        #coords1 = array([tuple(self.a['C3\''].get_vector()), tuple(self.b['C3\''].get_vector())], 'f')
        #coords2 = array([tuple(other.a['C3\''].get_vector()), tuple(other.b['C3\''].get_vector())], 'f')

        diff = coords1 - coords2 
        l = coords1.shape[0]
        r1 = sqrt(sum(sum(diff*diff))/l)

        return r1

    def reset_coord(self):
        """hack"""
        self.a = self.a_bk.copy()
        self.b = self.b_bk.copy()
        del self.a_bk
        del self.b_bk
        self.atoms_ab = self.get_atoms_ab()
        self.atoms_ba = self.get_atoms_ba()


if '__main__' == __name__:
    parser_pdb = PDB.PDBParser()
    struct = parser_pdb.get_structure('', 'test_data/R_s1t1r7.pdb')
    model = struct[0]
    chain = model['A']

    # 1 43
    a = chain[3]
    b = chain[43]
    print((a,b))


    bp = BasePair(a,b,struct)
    print(('coord:', bp.calc_coord()))
    print((bp.name))
    #dist = bg - bg2

    a = chain[2]
    b = chain[42]
    print((a,b))
    bp2 = BasePair(a,b, struct)
    print(('coord:', bp2.calc_coord()))
    print((bp2.name))

    print(('bp, bp2, distance', bp, bp2, bp2 - bp))

    print(('rmsd:', bp.calc_rmsd_to(bp2)))
    bp.overlap_with(bp2)
    bp.save()
    print(('overlapped rmsd', bp.calc_rmsd_to(bp2)))
    bp.reset_coord()
    print(('reset & rmsd', bp.calc_rmsd_to(bp2)))

    #print bp.get_atoms()
