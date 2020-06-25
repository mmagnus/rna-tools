
from rna_tools.tools.mini_moderna3.moderna.util.Errors import ModernaResidueError
from rna_tools.tools.mini_moderna3.moderna.Constants import BACKBONE_RIBOSE_ATOMS
from Bio.PDB.vectors import rotaxis2m
import math

def rotate_chi(resi, angle=90):
        """
        Rotates the base around the glycosidic bond
        by the given chi angle
        """
        try:
            vC1 = resi["C1'"].get_vector()
            vN = resi["N*"].get_vector()
        except KeyError:
            raise ModernaResidueError("Residue %d: could not find atom C1' or N9 (if A or G)or N1 (if C or U)." % resi.number )

        v_glycosidic = vN - vC1
        matrix = rotaxis2m(math.radians(angle), v_glycosidic)

        for atom in resi:
            if atom.name not in BACKBONE_RIBOSE_ATOMS and atom != resi["N*"]:
                v_atom = resi[atom.name].get_vector()
                vec = v_atom - vC1
                vec = vec.left_multiply(matrix)
                vec += vC1
                atom.set_coord([vec[0], vec[1], vec[2]])
