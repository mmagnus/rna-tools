#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .PDBFile import PDBFile


class RosettaPDBFile(PDBFile):
    """Class for preparing PDB files for ROSETTA/FARNA
    """
    def make_rna_rosetta_ready(self):
        """ROSETTA has some special requirements for PDB files.
    
        Arguments:
          * pdb_string = contents of PDB file as a string
    
        Output:
          * new PDB returned as a string
        """
        #self.check_and_get_first_model()
        self.resname_check_and_3to1()
        #self.terminate_chains()
        #self.remove_short_chains()
        result = []
        for l in self.pdb_lines:
            ## fix G__ -> __G
            if l.startswith('ATOM') and l[17] in ('A', 'U', 'C', 'G'):
                l = l[:16] + '   ' + l[17] + ' ' + l[21:]
            ##

            if l.startswith('ATOM') and l[17] in ('A', 'U', 'C', 'G'):
                line = l.replace('\'', '*')
                result.append(line[:18] + 'r' + line[19:])

            if l.startswith('ATOM') and l[19] in ('A', 'U', 'C', 'G'):
                line = l.replace('\'', '*')
                result.append(line[:18] + 'r' + line[19:])
            else:
                if not l.startswith('HETATM'):
                    result.append(l)
        self._apply_fix('make_rna_rosetta_ready', '\n'.join(result), result)
    
    
    def make_rna_non_rosetta(self):
        """Reverses changes required for ROSETTA, especially X -> rX.
    
        Arguments:
          * pdb_string = contents of PDB file as a string
    
        Output:
          * new PDB returned as a string
        """
        result = []
        for l in self.pdb_lines:
            if l.startswith('ATOM') and l[19] in ('A', 'U', 'C', 'G'):
                line = l.replace('\'', '*')
                result.append(line[:18] + ' ' + line[19:])
            else:
                result.append(l)
        self._apply_fix('make_rna_non_rosetta_ready', '\n'.join(result), result)