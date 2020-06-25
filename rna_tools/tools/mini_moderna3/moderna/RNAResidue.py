"""
Wrapper for PDB.Residue objects.
"""

import re
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from numpy import array
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaAlphabet import alphabet
from rna_tools.tools.mini_moderna3.moderna.analyze.BaseRecognizer import BaseRecognizer, BaseRecognitionError
from rna_tools.tools.mini_moderna3.moderna.analyze.BasePairCalculator import base_pair_calc
from rna_tools.tools.mini_moderna3.moderna.builder.CoordBuilder import build_coord

from rna_tools.tools.mini_moderna3.moderna.util.Errors import RNAResidueError

from rna_tools.tools.mini_moderna3.moderna.Constants import STANDARD_BASES, ANY_RESIDUE
from rna_tools.tools.mini_moderna3.moderna.Constants import PURINE_NEIGHBOR_TABLE, PYRIMIDINE_NEIGHBOR_TABLE
from rna_tools.tools.mini_moderna3.moderna.Constants import DONORS, ACCEPTORS
from rna_tools.tools.mini_moderna3.moderna.Constants import H_GENERATE_STEP, H_COVALENT_BOND
from rna_tools.tools.mini_moderna3.moderna.Constants import H_ANGLE_ONE, H_ANGLE_TWO


H_GENERATE_TORSIONS = list(range(0, 360, H_GENERATE_STEP))


class RNAResidue(Residue):
    """
    Supplements Bio.PDB.Residue object with functions to
    manage RNA-specific features.
    """
    br = BaseRecognizer()

    def __init__(self, pdb_residue, alphabet_entry=None, new_atoms=True):
        """
        Arguments:
        - residue as a Bio.PDB.Residue instance
        - optional: AlphabetEntry instance, if it is not given, the residue
                    is identified using BaseRecognizer (slow)
        """
        Residue.__init__(self, pdb_residue.id, pdb_residue.resname, '    ')
        self.number = pdb_residue.id[1]
        self.disordered = pdb_residue.disordered
        self.has_double_coord = self.__check_double_coord(pdb_residue)
        self.modified = None
        self.long_abbrev = None
        if alphabet_entry:
            abbrev = alphabet_entry.long_abbrev
        else:
            try:
                abbrev = self.br.identify_resi(pdb_residue)
            except BaseRecognitionError:
                abbrev = alphabet.get_short_original(ANY_RESIDUE)
        self.change_name(abbrev)

        self.identifier = str(self.id[1]).strip() + self.id[2].strip()
        self.__create_atoms(pdb_residue, new_atoms)

        # caches for H-bond calculation --> faster
        self._donors = None
        self._acceptors = None
        self._donor_hydrogens = {}

    def __check_double_coord(self, resi):
        """
        Checks whether any atoms in residue
        have alternative coordinates given in the pdb file.
        """
        if not self.disordered:
            return False
        for atom in resi:
            if atom.is_disordered():
                if len(atom.disordered_get_list()) > 1:
                    return True
        return False

    def __create_atoms(self, pdb_residue, new_atoms):
        if new_atoms:
            # copy all atoms, in case the original is manipulated.
            for atom in pdb_residue.child_list:
                if not atom.name[0] in '*H123':  # .startswith('H'):
                    element = re.sub('[\s\d]', '', atom.name)[0] or 'C'
                    new_at = Atom(atom.name, atom.coord, atom.bfactor,
                                  atom.occupancy, atom.altloc, atom.fullname,
                                  atom.serial_number, element=element)
                    self.add(new_at)
        else:
            # use the old atoms (saves time)
            [self.add(atom) for atom in pdb_residue.child_list]

    def __len__(self):
        """Returns number of atoms."""
        return len(self.child_list)

    def __repr__(self):
        """Returns string representation"""
        return '<Residue %s %s>' % (self.identifier, self.long_abbrev)

    def __getitem__(self, name):
        """
        Returns an atom like PDB.Residue,
        but interprets N* as the glycosidic N.
        """
        # KR: added to take care of the N1/N9 locally (important for LIR)
        # C,U ---> N1;    A,G ---> N9;    pseudouridine ---> C5;
        if name == 'N*':
            if self.long_abbrev in ['Y', 'm1acp3Y', 'Ym', 'm3Y']:
                return Residue.__getitem__(self, 'C5')
            if self.pyrimidine:
                return Residue.__getitem__(self, 'N1')
            elif self.purine:
                return Residue.__getitem__(self, 'N9')
            elif self.original_base == 'X':
                if 'N9' in self.child_dict:
                    return Residue.__getitem__(self, 'N9')
                elif 'N1' in self.child_dict:
                    return Residue.__getitem__(self, 'N1')
                else:
                    raise RNAResidueError('Cannot decide which atom to use for glycosidic N in residue %s' % self)
            elif 'N1' in self.child_dict:
                return Residue.__getitem__(self, 'N1')
            elif 'N9' in self.child_dict:
                return Residue.__getitem__(self, 'N9')
            else:
                raise RNAResidueError('Cannot decide which atom to use for glycosidic N in residue %s' % self)
        else:
            return Residue.__getitem__(self, name)

    def change_number(self, new_number):
        """Changes a residues number to the given string."""
        try:
            num = int(new_number.strip())
            self.id = (self.id[0], num, ' ')
        except:
            try:
                letter = new_number.strip()[-1]
                num = int(new_number.strip()[:-1])
                self.id = (self.id[0], num, letter)
            except ValueError:
                raise RNAResidueError('Invalid residue number: %s' % new_number)
        self.number = num
        self.identifier = new_number.strip()

    def change_name(self, new_name):
        """
        Changes the residues name.
        to a new name (as a long abbreviation if modified)
        """
        if new_name not in alphabet:
            new_name = alphabet.get_short_original(ANY_RESIDUE).long_abbrev
        aentry = alphabet[new_name]
        self.resname = aentry.pdb_abbrev
        self.long_abbrev = aentry.long_abbrev

        if new_name in STANDARD_BASES:
            # standard base
            self.modified = False
            self.id = (' ', self.id[1], self.id[2])
        elif aentry.original_base.upper() in "X":
            # unknown residues --> water, ions.
            self.modified = False
        else:
            # modified base
            self.modified = True
            self.id = ('H_' + aentry.pdb_abbrev, self.id[1], self.id[2])
            if aentry.pdb_abbrev == 'UNK':
                abbrev = '0' * (3 - len(aentry.new_abbrev)) + aentry.new_abbrev
                self.resname = abbrev
                self.id = ('H_' + abbrev, self.id[1], self.id[2])

        self._clear_caches()

    def _clear_caches(self):
        """Delete internally saved shortcuts"""
        self._donors = None
        self._acceptors = None
        self._donor_hydrogens = {}

    @property
    def alphabet_entry(self):
        """Returns an alphabet entry for this residue."""
        return alphabet[self.long_abbrev]

    @property
    def purine(self):
        """Returns True if the residue is a purine."""
        if self.original_base in ("G", "A"):
            return True

    @property
    def pyrimidine(self):
        """Returns True if the residue is a pyrimidine."""
        if self.original_base in ("C", "U"):
            return True

    @property
    def original_base(self):
        """Returns the unmodified base abbreviation."""
        return self.alphabet_entry.original_base

    @property
    def short_abbrev(self):
        """Returns a one-letter abbreviation of the residue."""
        return self.alphabet_entry.short_abbrev

    @property
    def new_abbrev(self):
        """Returns the Modomics nomenclature abbreviation."""
        return self.alphabet_entry.new_abbrev

    @property
    def pdb_abbrev(self):
        """Returns a three-letter PDB abbreviation."""
        return self.alphabet_entry.pdb_abbrev

    @property
    def full_name(self):
        """Returns the full name of the nucleotide."""
        return self.alphabet_entry.full_name

    @property
    def category(self):
        """Returns the cathegory of the nucleotide."""
        return self.alphabet_entry.category

    # ------------------- helper functions to work with atoms ---------------
    def get_atom_vector(self, name):
        """returns a vector for the given atom. N* encodes the glyco-N"""
        try:
            return self[name].get_vector()
        except KeyError:
            raise RNAResidueError('There is no atom %s in residue %s'
                                  % (name, self.identifier))

    def get_atoms_by_names(self, names, strict=False):
        """Generates atoms from the given list of names."""
        for name in names:
            try:
                yield self[name]
            except KeyError:
                if strict:
                    raise KeyError("Atom %s not found" % name)

    def check_atoms(self, names):
        """Returns True if all atom names in the given list exist."""
        try:
            [atom for atom in self.get_atoms_by_names(names, strict=True)]
            return True
        except KeyError:
            return False

    def get_bp(self, resi2):
        """returns an interaction type between two residues
        or None if there is no interaction"""
        return(base_pair_calc(self, resi2))

    # ------------ helper methods for h-bond calculation ----------------------
    def get_hbond_donors(self):
        """Generates atoms that are H-bond donors of this residue."""
        if self._donors:
            return self._donors
        key = self.original_base.strip()
        result = list(self.get_atoms_by_names(DONORS.get(key, [])))
        self._donors = result
        return result

    def get_hbond_acceptors(self):
        """Generates atoms that are H-bond acceptors of this residue."""
        if self._acceptors:
            return self._acceptors
        key = self.original_base.strip()
        result = list(self.get_atoms_by_names(ACCEPTORS.get(key, [])))
        self._acceptors = result
        return result

    def get_neighbors(self, atom):
        """Returns a list of atoms in the same residue connected by bonds."""
        result = []
        if 'N9' in self.child_dict:
            nb_dict = PURINE_NEIGHBOR_TABLE
        else:
            nb_dict = PYRIMIDINE_NEIGHBOR_TABLE
        for name in nb_dict.get(atom.fullname):
            child = self.child_dict.get(name)
            if child:
                result.append(child)
        return result

    def get_donor_hydrogens(self, donor):
        """
        Returns a list of coord records of hypothetical hydrogen positions.
        If the donor has two neighbors, this will be a single position, if
        it has only one, a rotation will be performed in 10 degree steps.
        Atoms with 3 or more neighbors will be rejected.
        """
        # TODO: refactor this out.
        if donor.name in self._donor_hydrogens:
            return self._donor_hydrogens[donor.name]

        hydrogens = []
        neighbors = self.get_neighbors(donor)
        don_vec = donor.get_vector()
        sup_vec1 = None  # coordinates next to donor
        sup_vec2 = None  # coordinates to calc torsion

        if len(neighbors) == 1:
            sup_vec1 = neighbors[0].get_vector()
            neighbors2 = self.get_neighbors(neighbors[0])
            sup_vec2 = None
            while neighbors2 and sup_vec2 is None:
                next = neighbors2.pop()
                if next != donor:
                    sup_vec2 = next.get_vector()
            # bad case: no neighbors to construct 2nd support vec
            if sup_vec2 is None:
                sup_vec2 = (don_vec ** sup_vec1)
            angle = H_ANGLE_ONE
            torsions = H_GENERATE_TORSIONS

        elif len(neighbors) == 2:
            sup_vec1 = neighbors[0].get_vector()
            sup_vec2 = neighbors[1].get_vector()
            angle = H_ANGLE_TWO
            torsions = [180.0]

        if sup_vec1 is not None and sup_vec2 is not None:
            # create hydrogen positions
            for torsion in torsions:
                h_pos = build_coord(sup_vec2, sup_vec1, don_vec,
                                    H_COVALENT_BOND, angle, torsion)
                h_pos = array([h_pos[0], h_pos[1], h_pos[2]])
                hydrogens.append(h_pos)

        self._donor_hydrogens[donor.name] = hydrogens
        # self.write_hydrogens(hydrogens)
        return hydrogens
        # TODO: what if there are 0 neighbors (water)?
