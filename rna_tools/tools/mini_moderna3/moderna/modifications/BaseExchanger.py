
from rna_tools.tools.mini_moderna3.moderna.modifications.ResidueEditor import ResidueEditor
from rna_tools.tools.mini_moderna3.moderna.modifications.ModificationRemover import remove_modification
from rna_tools.tools.mini_moderna3.moderna.util.Errors import ExchangeBaseError
from rna_tools.tools.mini_moderna3. moderna.util.LogFile import log
from rna_tools.tools.mini_moderna3.moderna.Constants import UNKNOWN_RESIDUE_SHORT, B_FACTOR_EXCHANGE, BASE_PATH


class BaseExchanger(ResidueEditor):

    def get_exchange_rule(self, resi,new_name):
        """
        Prepares a rule for exchanging a base.
        Returns a dict with the rule.
        {
        'fixed':[link atom names from original residue],
        'moved':[link atom names from new base],
        'remove':[atom names that must be removed from the original residue (old base atoms)]
        }
        """
        rule = {}
        if resi.purine:
            rule['fixed'] = ['N9', 'C4', 'C8']
            rule['remove'] = ['N9', 'C8', 'C4', 'N1', 'C6', 'C5', 'N7', 'C8', 'N3', 'O6', 'N2', 'N6', 'C2']
        elif resi.pyrimidine:
            rule['fixed'] = ['N1', 'C2', 'C6']
            rule['remove'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'N4', 'C5', 'C6']
        else:
            raise ExchangeBaseError('Residue %s: could not get exchange rule for name %s' % resi.identifier, new_name)

        if new_name in ['A', 'G']:
            rule['moved'] = ['N9', 'C4', 'C8']
        elif new_name in ['C', 'U']:
            rule['moved'] = ['N1', 'C2', 'C6']
        else:
            raise ExchangeBaseError('Residue %s: could not get exchange rule for name %s' % resi.identifier, new_name)
        return rule


    def exchange_base(self, resi, new_name):
        """
        Exchanges standard bases in a residue.

        Arguments:
        - a new base name ('A','G','C' or 'U')
        """
        if resi.long_abbrev == UNKNOWN_RESIDUE_SHORT:
            self.mutate_unknown_residue()
        if resi.modified:
            remove_modification(resi)

        rule = self.get_exchange_rule(resi, new_name)
        new_base = self.parse.get_structure(resi.original_base, BASE_PATH+new_name+'.ent')[0]['C'][(' ', 54, ' ')]

        self.superimpose.get_atoms([resi], rule['fixed'], 'fixed')
        self.superimpose.get_atoms([new_base], rule['moved'], 'moved')
        self.superimpose.moved_atoms = new_base.child_list
        self.superimpose.superimpose()

        for atom in rule['remove']:
            if atom in resi.child_dict.keys():
                resi.detach_child(atom)
        for atom in new_base:
            resi.add(atom)
        resi.change_name(new_name)
        self.set_bfactor(resi, B_FACTOR_EXCHANGE)



def exchange_base(resi, new_name):
    """
    Exchanges base in given residue.

    Arguments:
    - residue
    - new residue name (A, G, C or U)
    """
    old_name = resi.long_abbrev
    bex = BaseExchanger()
    bex.exchange_base(resi, new_name)
    log.write_message('Residue %s: base exchanged (%s ---> %s), residue added to model.' %(resi.identifier, old_name, new_name))



