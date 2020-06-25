#!/usr/bin/env python
#
# StructureLibrary.py
#
# Cache for ModernaStructure obiects.
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

from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.ModernaFragment import ModernaFragment
import os, tempfile

class StructureLibraryError(Exception): pass

class StructureLibrary:
    """
    A class that produces ModernaStructure objects.
    If the same PDB file is loaded twice, the structure is cached.
    The usage of StructureLibrary requires that the PDB files do not change 
    on disk at runtime.
    (Implements the Factory Pattern for ModernaStructure)
    """
    def __init__(self, path):
        self.path = path
        self.structures = {}
        self.sequences = {}

    def get_structure(self, name, chain='A'):
        """Returns a ModernaStructure object from a PDB file."""
        key = (name, chain)
        if self.structures.has_key(key):
            seq = self.sequences[key]
            struc = ModernaStructure('file', self.path+name, chain, seq, new_atoms=False)
            # KR: the following lines are SLOWER
            #struc = self.structures[key]
            #struc = ModernaStructure('residues', [r for r in struc], seq)
        else:
            seq = self.sequences.get(key, None)
            struc = ModernaStructure('file', self.path+name, chain, seq, new_atoms=False)
            seq = struc.get_sequence()
            
            self.structures[key] = struc
            self.sequences[key] = seq
            struc = ModernaStructure('residues', [r for r in self.structures[key]], seq)
        return struc
        
    def find_resi_in_lines(self, lines, resi):
        """Performs binary tree search for the residue and returns first and last position."""
        # convert residue number to number/letter pair
        if resi[-1] in '0123456789':
            res_int = int(resi)
            res_char = ' '
        else:
            res_int = int(resi[:-1])
            res_char = resi[-1]
        # perform the tree search
        found = False
        step = len(lines)/2
        i = step
        while not found and step>0:
            step = step/2
            # compare residue numbers
            probe = int(lines[i][22:26].strip())
            if probe == res_int: found = True
            elif int(probe)<res_int: i += step
            else: i -= step
        if not found:
            raise StructureLibraryError("Could not find residue '%s'"%resi)
            
        # found a position, now determine begin and end
        begin, end = i, i
        while begin>0 and int(lines[begin-1][22:26].strip()) == res_int:
            begin -= 1
        while end<len(lines) and int(lines[end][22:26].strip()) == res_int:
            end += 1
        return begin, end
        
    def create_pdb_sniplet_file(self, name, chain, resi_from, resi_to, sequence=None):
        """Returns filename of excerpt from PDB file containting only certain residues."""
        lines = open(self.path+name)
        chainlines = [line for line in lines if len(line)>21 and line[21]==chain] 
        chainlines = [line for line in chainlines if line.startswith('ATOM') or line.startswith('HETATM')]
        start, none = self.find_resi_in_lines(chainlines, resi_from)
        none, end = self.find_resi_in_lines(chainlines, resi_to)
        out = chainlines[start:end]
        tmpname = tempfile.mktemp()
        open(tmpname, 'w').writelines(out)
        return tmpname

    def get_structure_part(self, name, chain, resi_from, resi_to, sequence=None):
        """Reads only certain residues from a PDB file to a structure object."""
        tmpname = self.create_pdb_sniplet_file(name, chain, resi_from, resi_to, sequence)
        struc = ModernaStructure('file', tmpname, chain, seq=sequence)
        os.remove(tmpname)
        return struc
        
    def get_fragment_part(self, name, chain, resi_from, resi_to, \
                        anchor5, anchor3, sequence, keep, seq):
        """Returns a ModernaFragment53 from the library."""
        struc = self.get_structure_part(name, chain, resi_from, resi_to, seq)
        return ModernaFragment53(struc, anchor5, anchor3, sequence, keep=keep)
     

library = StructureLibrary('')
