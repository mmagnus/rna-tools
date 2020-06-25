
# analyzes whether the backbones of RNA residues are connected.

from rna_tools.tools.mini_moderna3.moderna.analyze.GeometryParameters import BACKBONE_DIST_MATRIX, \
    PHOSPHATE_DIST_MATRIX, O3_P_DIST_HI
from rna_tools.tools.mini_moderna3.moderna.Constants import BACKBONE_ATOMS,  \
    BACKBONE_RIBOSE_ATOMS_WITHOUT_O2

DIST_TOLERANCE = 1.05

# distance for intra-residue backbone clashes
CONGESTION_DISTANCE = 1.5

ATOMS_5P = ["P", "O5'", "C5'", "C4'"]
ATOMS_3P = ["C3'", "O3'"]

BB_SET_5P = ["P", "OP1", "OP2", "C5'", "O5'"]
BB_SET_3P = ["O3'"]
BB_SET = BB_SET_5P + BB_SET_3P

BONDS = {
        "P":["OP1", "OP2", "O5'"],
        "OP1":["P"],
        "OP2":["P"],
        "C5'":["C4'", "O5'"],
        "O3'":["C3'"],
        "O5'":["C5'", "P"],
        }


def are_residues_connected(res1, res2):
    """
    Checks whether two residues are connected.
    Distances on the backbone are within norm + tolerance.

    Arguments:
    * upstream residue as RNAResidue object
    * downstream residue as RNAResidue object
    """
    try:
        if res1["O3'"] - res2["P"] > O3_P_DIST_HI * DIST_TOLERANCE:
            return False
        for tup in [("C3'", "O3'"), ("C4'", "C3'")]:
            a1, a2 = tup
            low, hi = BACKBONE_DIST_MATRIX[tup]
            if res1[a1] - res1[a2] > hi * DIST_TOLERANCE:
                return False
        for tup in [   ("P", "O5'"), ("O5'", "C5'"), ("C5'", "C4'")]:
            a1, a2 = tup
            low, hi = BACKBONE_DIST_MATRIX[tup]
            if res2[a1] - res2[a2] > hi * DIST_TOLERANCE:
                return False
        return True
    except KeyError:
        # missing atoms
        return False


def is_chain_continuous(chain):
    """
    Checks whether a chain is continuous.
    Check whether all subsequent pairs of residues are connected.

    Arguments:
        chain - RNAChain object
    """
    keys_sorted = []
    for resi in chain:
        keys_sorted.append(resi.identifier)

    pairs_to_check = [(keys_sorted[x], keys_sorted[x+1]) for x in range(len(keys_sorted)-1)]

    for pair in pairs_to_check:
        if not are_residues_connected(chain[pair[0]], chain[pair[1]]):
            return False
    return True


# ------------- METHODS FOR CHECKING RESIDUE INTEGRITY --------------------
def is_backbone_complete(self):
    """Returns True if all backbone atoms are present."""
    return self.check_atoms(BACKBONE_ATOMS)


def is_ribose_complete(self):
    """Returns True if all ribose atoms are present."""
    return self.check_atoms(BACKBONE_RIBOSE_ATOMS_WITHOUT_O2)


def is_backbone_intact(self, tolerance=1.0, mode=None):
    """
    Checks whether all the backbone atoms in the residue are connected.
    Returns True/False value. In case any backbone atoms are missing\
    raises an exception.
    """
    if mode == "5'":
        atoms = ATOMS_5P
    elif mode == "3'":
        atoms = ATOMS_3P
    try:
        for atom in BACKBONE_DIST_MATRIX:
            if mode and (atom[0] not in atoms and atom[1] not in atoms):
                continue
            low_dist, hi_dist = BACKBONE_DIST_MATRIX[atom]
            dist = self[atom[0]] - self[atom[1]]
            if not (low_dist <= dist <= hi_dist * tolerance):
                return False
        return True
    except KeyError: # missing atoms
        return False


def is_phosphate_intact(self, tolerance=1.0):
    """Checks whether P-OP1 and P-OP2 distances are OK."""
    try:
        for atom in PHOSPHATE_DIST_MATRIX:
            low_dist, hi_dist = PHOSPHATE_DIST_MATRIX[atom]
            dist = self[atom[0]] - self[atom[1]]
            if not (low_dist <= dist <= hi_dist * tolerance):
                return False
        return True
    except KeyError: # missing atoms
        return False


def is_backbone_congested(self, congestion_dist=CONGESTION_DISTANCE, \
                    mode=None):
    """Checks whether backbone atoms clash into others."""
    atoms = BB_SET
    if mode == "5'":
        atoms = BB_SET_5P
    elif mode == "3'":
        atoms = BB_SET_3P
    for bb_name in atoms:
        try:
            atom1 = self[bb_name]
            for atom2 in self:
                if atom2.name in atoms:
                    continue
                if atom2.name in BONDS[bb_name]:
                    continue
                dist = atom2-atom1
                if dist < congestion_dist:
                    return True
        except KeyError:
            pass # skip missing atoms.
