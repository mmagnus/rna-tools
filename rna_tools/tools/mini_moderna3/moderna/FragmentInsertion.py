#!/usr/bin/env python
#
# FragmentInsertion
#
# Tool to insert fragments into structures.
# 
# http://iimcb.genesilico.pl/moderna/ 
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

class FragmentInserter(object):
    """
    Inserts fragments into models.
    """
    def prepare_fragment(self, fragment, model):
        """Prepares the fragment for insertion."""
        fragment.superimpose()
        fragment.prepare_anchor_residues() 
        fragment.renumber(model) 
        fragment.apply_seq() 

    def insert_fragment(self, fragment, model):
        """Inserts a ModernaFragment instance into the model."""
        self.prepare_fragment(fragment, model)
        fragment.fix_backbone()
        for resi in fragment.get_resi_to_remove(model): 
            model.remove_residue(resi)
        for resi in fragment.struc: 
            model.add_residue(resi)

