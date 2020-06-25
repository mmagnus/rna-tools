#!/usr/bin/env python
#
# Template.py
#
# Template structure for building a RNA model.
#
# http://iimcb.genesilico.pl/moderna/ 
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Magdalena Rother"
__email__ = "rother.magdalena@gmail.com"
__status__ = "Production"

from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure

class Template(ModernaStructure):
    """
    Wrapper for template structure.

    Arguments:
    * data: structure/chain/file name/residues list
    * data type specyfication ('structure'/'chain'/'file'/'residues')
    * chain name
    """
    def __init__(self, template_data, template_data_type='file', \
                 template_chain_name='A', seq=None):
        #TODO: same order of arguments as in ModernaStructure
        ModernaStructure.__init__(self, data_type= template_data_type, \
                            data=template_data, chain_name=template_chain_name, seq=seq)
        self.template_residues = {}
        self.set_template_numeration()

    def set_template_numeration(self):
        """ 
        Prepares dict of template residues. 
        The key for each template residue is the number of it residue in the sequence.
        The first residue in the template sequence has key '1' and so on.
        This numeration is independant from residues numeration in the template pdb structure.
        """        
        # connects ModernaResidues to indices in the alignment
        self.template_residues = {}
        number = 1
        for resi in self:
            self.template_residues[str(number)] = resi
            number += 1
        if self.check_letters_in_residue_numeration(): 
            self.renumber_chain()
