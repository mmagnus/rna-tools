"""
Removing modifications from RNA nucleotides
"""

from rna_tools.tools.mini_moderna3.moderna.modifications.ResidueEditor import ResidueEditor
from rna_tools.tools.mini_moderna3.moderna.util.Errors import RemoveModificationError
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log
from rna_tools.tools.mini_moderna3.moderna.Constants import BASE_PATH, BACKBONE_RIBOSE_ATOMS, B_FACTOR_REMOVE_MODIF


class ModificationRemover(ResidueEditor):

    def get_remove_rule(self, resi):
        """
        Prepares a rule for removing modification from a residue.
        Returns a dict with a rule:
        {
        'fixed': [link atom names from an original residue],
        'moved': [link atom names from a new (standard) base]
        }
        """
        if resi.long_abbrev in ['Y', 'Ym', 'm3Y', 'm1Y', 'm1acp3Y']:
            return {'fixed':['C5', 'C4', 'C6'], 'moved':['N1', 'C2', 'C6']}
        elif resi.purine:
            return {'fixed':['N9', 'C8', 'C4'], 'moved':['N9', 'C8', 'C4']}
        elif resi.pyrimidine:
            return {'fixed':['N1', 'C2', 'C6'], 'moved':['N1', 'C2', 'C6']}
        else:
            raise RemoveModificationError('Residue %s: could not get a removing rule.' % resi.identifier)


    def remove_modification(self, resi):
        """
        Removes a modification from this residue.
        It removes all unnecessary atoms and adds a standard base,
        corresponding to the originating base of the modified one.
        according to removing modification rule (prepared by get_remove_rule()).
        """
        if not resi.modified:
            raise RemoveModificationError('Residue %s: the residue does not have any modification. Could not remove modification.'  % resi.identifier)
        elif resi.long_abbrev in ['X', 'Xm']:
            raise RemoveModificationError('Residue %s: unidentified residue. Could not remove modification.' % resi.identifier)
        elif resi.long_abbrev in ['dA', 'dG', 'dC', 'dT']:
            rule = {
                'modification_name': resi.long_abbrev,
                'original_base': resi.original_base,
                'remove': ''
            }
            if resi.long_abbrev == 'dT':
                if resi.child_dict.has_key('C7'):
                    rule['remove'] = 'C7'
                elif resi.child_dict.has_key('C5M'):
                    rule['remove'] = 'C5M'
            rule['moved_link_atoms'] = ["C3'", "C2'", "C1'"]
            rule['fixed_link_atoms'] = ["C3'", "C2'", "C1'"]
            rule['fragment_file_name'] = 'C1pC2pC3p_O2p.pdb'
            rule['pdb_abbrev'] = 'D' + resi.long_abbrev[1]
            self.add_single_fragment(resi, rule)
        else:
            struc = self.parse.get_structure(resi.original_base, BASE_PATH + resi.original_base + '.ent')
            new_base = struc[0]['C'][(' ', 54, ' ')]
            triplet_names = self.get_remove_rule(resi)

            self.superimpose.get_atoms([resi], triplet_names['fixed'], 'fixed')
            self.superimpose.get_atoms([new_base], triplet_names['moved'], 'moved')
            self.superimpose.moved_atoms = new_base.child_list
            self.superimpose.superimpose()

            try:
                for atom in resi.child_list[:]:
                    if atom.id not in BACKBONE_RIBOSE_ATOMS:
                        resi.detach_child(atom.id)
                for atom in new_base:
                    resi.add(atom)
            except:
                raise RemoveModificationError('Residue %s: could not remove unnecessary and add proper atoms' % resi.identifier)

        resi.change_name(resi.original_base)
        resi.modified = False
        self.set_bfactor(resi, B_FACTOR_REMOVE_MODIF)


def remove_modification(resi):
    """Removes modification from a residue."""
    old_name = resi.long_abbrev
    m = ModificationRemover()
    m.remove_modification(resi)
    log.write_message('Residue %s: modification removed (%s ---> %s).' %(resi.id, old_name, resi.long_abbrev))


def remove_all_modifications(struc):
    """Removes all modifications from a structure."""
    for resi in struc:
        if resi.modified:
            remove_modification(resi)
