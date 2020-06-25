"""
Add modifications to a residue
"""

from rna_tools.tools.mini_moderna3.moderna.modifications.ResidueEditor import ResidueEditor
from rna_tools.tools.mini_moderna3.moderna.modifications.BaseExchanger import BaseExchanger
from rna_tools.tools.mini_moderna3.moderna.modifications.ModificationRemover import ModificationRemover
from rna_tools.tools.mini_moderna3.moderna.util.Errors import AddModificationError
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log
from rna_tools.tools.mini_moderna3.moderna.Constants import ANY_RESIDUE, MISSING_RESIDUE, \
    UNKNOWN_RESIDUE_SHORT, B_FACTOR_ADD_MODIF, \
    ADDING_MODIFICATION_RULES_PATH


def parse_modification_rules(separator=' | '):
    """
    Prepares a rule for adding a modification.
    Rules describe which fragments add and how to do this
    to obtain a residue with given modification.
    Returns dict of list of dicts with rules for adding a single fragment.

    Keys in each rule dict: ['modification_name', 'original_base', 'remove',
    'moved_link_atoms', 'fixed_link_atoms', 'fragment_file_name', 'pdb_abbrev']
    """
    rules = {}
    try:
        infile = open(ADDING_MODIFICATION_RULES_PATH)
    except IOError:
        log.write_message('File does not exist: %s ' % ADDING_MODIFICATION_RULES_PATH)
        return {}

    for line in infile:
        line = line.strip().split(separator)
        if len(line) >= 7:
            mod_name = line[0].strip()
            rules.setdefault(mod_name, [])
            rule = {}
            rule['modification_name'] = line[0]
            rule['original_base'] = line[1]
            rule['remove'] = line[2]
            rule['moved_link_atoms'] = line[3].split(',')
            rule['fixed_link_atoms'] = line[4].split(',')
            rule['fragment_file_name'] = line[5]
            rule['pdb_abbrev'] = line[6]
            rules[mod_name].append(rule)
    return rules

MODIFICATION_RULES = parse_modification_rules()


class ModificationAdder(ResidueEditor):

    def add_modification(self, resi, modification_name):
        """
        Adds a modification to a residue.
        It adds single fragments (add_single_fragment)
        according to adding modification rules (get_modification_rules).

        Arguments:
        - modification name (as a long abbreviation)
        """
        try:
            if modification_name in [ANY_RESIDUE, MISSING_RESIDUE]:
                raise AddModificationError('Residue %s: expected a modification name, instead got missing/any residue abbreviation "%s"'\
                                % (resi.identifier, modification_name))
            else:
                if resi.long_abbrev == UNKNOWN_RESIDUE_SHORT:
                    self.mutate_unknown_residue(resi)
                if resi.modified:
                    rem = ModificationRemover()
                    rem.remove_modification(resi)
                rules = MODIFICATION_RULES.get(modification_name, [])
                if not rules:
                    raise AddModificationError('Residue %s: there is no rule for adding this modification. Check modification name "%s".' \
                        %(resi.identifier, modification_name))
                else:
                    if rules[0]['original_base'] != resi.original_base:
                        bex = BaseExchanger()
                        bex.exchange_base(resi, rules[0]['original_base'])
                    for rule in rules:
                        self.add_single_fragment(resi, rule)
                    resi.change_name(modification_name)
                    self.set_bfactor(resi, B_FACTOR_ADD_MODIF)
        except IOError:
            raise AddModificationError('Residue %s: could not add modification.' % resi.identifier)


def add_modification(resi, long_abbrev):
    """Adds modification with given abbreviation"""
    old_name = resi.long_abbrev
    add = ModificationAdder()
    add.add_modification(resi, long_abbrev)
    log.write_message('Residue %s: modification added (%s ---> %s).' %(resi.identifier, old_name, long_abbrev))


