#!/usr/bin/env python
#
# IsostericityModelling.py
#
# Module for modeling isosteric base pairs.
#
# http://iimcb.genesilico.pl/moderna/
#

__author__ = "Pawel Skiba, Magdalena Rother, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__maintainer__ = "Pawel Skiba"
__email__ = "pw.skiba@gmail.com"
__status__ = "Prototype"


#from rna_tools.tools.mini_moderna3.moderna.import *
from rna_tools.tools.mini_moderna3.moderna.ModernaResidue import ModernaResidue
from IsostericityMatrices import IsostericityMatrices
from Isostericity import Isostericity
from rna_tools.tools.mini_moderna3.moderna.util.Errors import IsostericityError
import os

def create_alignment(model, template):
    """
    Creates output file with sequences in fasta format
    template and model are tuples containing pdb file name and chain id
    """
    file_name = "align_"+model[0][:-4]+"_"+template[0][:-4]+".fasta"
    alignment = open(file_name,"w")
    model_seq = str(ModernaStructure('file', model[0], model[1]).get_sequence())
    template_seq = str(ModernaStructure('file', template[0], template[1]).get_sequence())
    alignment.write(">"+model[0][:-4]+"\n"+model_seq+"\n")    
    alignment.write(">"+template[0][:-4]+"\n"+template_seq)    
    alignment.close()
    try:
        os.system("clustalw " + file_name + " -output=fasta")
    except:
        raise IsostericityError("clustalw must be installed to create an alignment")
    return file_name


class IsostericModelling:
    """
    Uses ModeRNA to create standard model and then models 
    isosteric base pairs using Isostericity.py
    """
    def __init__(self, template_file, template_chain, alignment):
        try:
            self.template = load_template(template_file, template_chain)
            self.alignment = load_alignment(alignment)
        except:
            raise IsostericityError("Wrong arguments")
        self.matrices = IsostericityMatrices()
        self.model = None
        self.modelling()


    def modelling(self):
        """
        Creates standard model and then starts isosteric bases modelling
        """
        self.model = create_model(self.template, self.alignment)        
        self.model.write_pdb_file("model_standard.pdb")               

        model_bp_id = self.model['79'].resname + self.model['98'].resname
        template_bp = (self.template['79'], self.template['98'])
        iso = Isostericity(template_bp, model_bp_id, 2.0, "trans", True)
        new_bp = iso.result_bp

        if not new_bp:
            raise IsostericityError("Base pair is not isosteric - see log for details")

        self.replace_bp_in_model(new_bp)    
        #self.model.fix_backbone()
        self.model.write_pdb_file("model_isosteric.pdb")


    def replace_bp_in_model(self, new_bp):
        """
        Replaces base pair in model using isosteric one
        """
        self.replace_single_base(new_bp[0], new_bp[0].identifier)
        self.replace_single_base(new_bp[1], new_bp[1].identifier)
        

    def replace_single_base(self, new_base, old_base):
        """
        Replaces single base
        """
        print "Modelling residue: "+str(new_base)
        self.model.add_residue(new_base, old_base, False)   
             


############################################################

if __name__ == '__main__':
#    fragments = (('1j5e_kt23.pdb','A'),('1s72_kt7.pdb', '0'))   #automatic alignment using clustalw
#    align = create_alignment(*fragments)
    align = "align_1j5e_kt23_1s72_kt7.fasta"
    model_args = ('1s72_kt7.pdb', '0',align)
    modelling = IsostericModelling(*model_args)

