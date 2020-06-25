#!/usr/bin/env python
#
# IsostericityMatrices.py
#
# Module for modeling isosteric base pairs.
#
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Pawel Skiba, Magdalena Rother, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Pawel Skiba"
__email__ = "pw.skiba@gmail.com"
__status__ = "Prototype"


from rna_tools.tools.mini_moderna3.moderna.Constants import DATA_PATH
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log

class IsostericityMatrices:
    """
    Isostericity matrices implementation in Python
    """
    def __init__(self):
        self.source = open(DATA_PATH+"IsostericityMatrices.txt","r")
        self.matrices = self.import_matrices()
        self.source.close()
        
    def import_matrices(self):
        """
        Reads data from source txt files and prepare the dictionary
        """
        matrices = {}
        for line in self.source:
            if line.startswith("Type: "):
                master_key = line.split(" ")[1]
                type_key = line.split(" ")[2].strip()
                if master_key not in matrices.keys():
                    matrices[master_key] = {}
                matrices[master_key][type_key] = {}
            elif line.startswith("\n"):
                continue
            else:
                data = line.split(": ")
                if data[1] == '\n':
                    matrices[master_key][type_key][data[0]] = None
                else:
                    matrices[master_key][type_key][data[0]] = float(data[1].strip()) 
        return matrices 
        
    def check_isostericity(self, bp1, bp2, interact_type, max_value=1.0):
        """
        Returns True if basepair1 is isosteric to basepair2 when interaction type is interact_type
        """
        try:
            result = self.matrices[bp1][interact_type][bp2] 
            log.write_message(bp1+"->"+bp2+" ("+interact_type+")")
            return result <= max_value and result != None
        except:
            log.write_message("No information in IsostericityMatrices about: "+bp1+"->"+bp2+" ("+interact_type+")")
            return False
        
        
    def show_isosteric_bp(self, bp1, interact_type, max_value=1.0):
        """
        Returns a tuple with all base pairs isosteric to bp1 when interaction type is interact_type
        """
        pairs = self.matrices[bp1][interact_type]
        result = [bp for bp in pairs if pairs[bp] <= max_value and pairs[bp] != None]
        return tuple(result)
        
                
##########################################################################                
                
if __name__ == '__main__':
    test = IsostericityMatrices()
    print test.check_isostericity('AC','GU','cWW') #Is isosteric
    print test.check_isostericity('AA','GG','cWW') #Is None
    print test.check_isostericity('AC','UU','cWW') #Is not isosteric
    print test.show_isosteric_bp('AC','cWW')
