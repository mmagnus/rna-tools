
from rna_tools.tools.mini_moderna3.moderna.modifications.ResidueEditor import ResidueEditor, make_backbone_only_residue
from rna_tools.tools.mini_moderna3.moderna.modifications.ModificationRemover import remove_modification
from rna_tools.tools.mini_moderna3.moderna.modifications.ModificationAdder import add_modification
from rna_tools.tools.mini_moderna3.moderna.modifications.BaseExchanger import exchange_base
from rna_tools.tools.mini_moderna3.moderna.Constants import UNKNOWN_RESIDUE_SHORT


class ResidueModifier(ResidueEditor):
    """
    Convenience class: change every residue into every other.
    Decides by itself which operation to perform.
    """

    def mutate(self, resi, abbrev):
        """Exchange base, add or removing modification."""
        if abbrev == UNKNOWN_RESIDUE_SHORT:
            make_backbone_only_residue(resi)
        if resi.long_abbrev != abbrev:
            if resi.long_abbrev == UNKNOWN_RESIDUE_SHORT:
                self.mutate_unknown_residue(resi)
            if resi.modified:
                remove_modification(resi)
            if abbrev not in ['A', 'U', 'G', 'C']:
                add_modification(resi, abbrev)
            else:
                exchange_base(resi, abbrev)


def modify_residue(resi, abbrev):
    """
    Convenience function exchanging base, adding and removing modification.

    Arguments:
    - a long abreviation ('A', 'mA')
    """
    mod = ResidueModifier()
    mod.mutate(resi, abbrev)
