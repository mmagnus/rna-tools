
from Bio.PDB import PDBParser
from rna_tools.tools.mini_moderna3.moderna.ModernaSuperimposer import ModernaSuperimposer
from rna_tools.tools.mini_moderna3.moderna.util.Errors import ModernaResidueError
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log
from rna_tools.tools.mini_moderna3.moderna.Constants import MODIFICATION_FRAGMENTS_PATH, \
    BACKBONE_RIBOSE_ATOMS, ANY_RESIDUE, RIBOSE


class ResidueEditor:
    """
    Deals with residue objects.
    Supplements RNAResidue object with functions to
    - exchange bases
    - add modifications
    - remove modifications
    - rotat chi angle
    The type of each ModernaResidue is automatically recognized, and has a couple
    of long and short names as attributes:
    - long_abbrev
    - short_abbrev (one letter abbreviation)
    - original_base
    - full_name
    - modified (True/False)
    """
    parse = PDBParser()
    superimpose = ModernaSuperimposer()


    def add_single_fragment(self, resi, rule):
        """
        Adds a fragment to a residue.

        Arguments:
        - an adding rule dict (prepared by get_modification_rules())
        """
        try:
            fragment = self.parse.get_structure('fragment', MODIFICATION_FRAGMENTS_PATH + rule['fragment_file_name'])[0]['A'][('H_UNK', 1, ' ')]
        except IOError:
            raise ModernaResidueError('File does not exist: %s' % MODIFICATION_FRAGMENTS_PATH)       

        self.superimpose.get_atoms([resi], rule['fixed_link_atoms'], 'fixed')
        self.superimpose.get_atoms([fragment], rule['moved_link_atoms'], 'moved')
        self.superimpose.moved_atoms = fragment.child_list
        self.superimpose.superimpose()

        if rule['remove']: 
            delete_atoms = rule['remove'].split(',') + rule['moved_link_atoms']
        else:
            delete_atoms = rule['moved_link_atoms']
        try:
            for atom_name in delete_atoms:
                resi.detach_child(atom_name)
            for atom in fragment:
                resi.add(atom)
        except ModernaResidueError as e:
            raise e
        except:
            raise ModernaResidueError('Residue %s: could not remove unnecessary and add proper atoms' % resi.identifier)
            #TODO: remove except:


    def mutate_unknown_residue(self, resi):
        """
        Makes a Cytosine out of unknown residue (X, .) on ribose and N* (N9,N1) atoms.
        When ribose and N* are nor present raises an error.
        C can be then changed into any modification / standard residues.
        """
        for x in RIBOSE:
            if not resi.has_id(x):
                raise ModernaResidueError('Residue %s: cannot mutate unknown residue. Missing ribose atom %s' %(resi.identifier, x))
        if not resi.has_id('N9') and not resi.has_id('N1'):
            raise ModernaResidueError('Residue %s: cannot mutate unknown residue. Missing ribose atom N* (N1 or N9)' %(resi.identifier))

        try:
            fragment = self.parse.get_structure('fragment', MODIFICATION_FRAGMENTS_PATH +'ribose_C.pdb')[0]['A'][(' ', 1, ' ')]
        except IOError:
            raise ModernaResidueError('File does not exist: %s' % MODIFICATION_FRAGMENTS_PATH+'ribose_C.pdb')

        self.superimpose.get_atoms([resi], ["O4'", "C1'", "C2'"], 'fixed')
        self.superimpose.fixed.append(resi["N*"])
        self.superimpose.get_atoms([fragment], ["O4'", "C1'", "C2'", 'N1'], 'moved')
        self.superimpose.moved_atoms = fragment.child_list
        self.superimpose.superimpose()

        for x in list(resi): # must copy items before deleting.
            if  x.name not in BACKBONE_RIBOSE_ATOMS:
               resi.detach_child(x.name)

        for x in fragment:
            if x.name not in  ["O4'", "C1'", "C2'"]:
                resi.add(x)

        resi.change_name('C')
        log.write_message('Residue %s: patched to construct base atoms from scratch.' %(resi.identifier))


    def set_bfactor(self, resi, b_value):
        """Sets the same b factor for all atoms in the residue"""
        for atom in resi:
            atom.set_bfactor(b_value)



def make_backbone_only_residue(resi, include_N=True):
    """
    Cuts the base out of a residue so only backbone and ribose atoms stay.
    """
    backbone_atoms = BACKBONE_RIBOSE_ATOMS[:]
    if include_N:
        nglyc = resi["N*"]
        backbone_atoms += [nglyc.name.strip()]
    for atom in resi.child_list[:]:
        # copying list because list is modified while iterated
        if atom.name.strip() not in backbone_atoms:
            resi.detach_child(atom.id)

    for at_name in backbone_atoms:
        if at_name not in list(resi.child_dict.keys()):
            raise ModernaResidueError('Residue %s: backbone is not complete. Missing atom %s' %(resi.identifier, at_name))
    resi.change_name(ANY_RESIDUE)
