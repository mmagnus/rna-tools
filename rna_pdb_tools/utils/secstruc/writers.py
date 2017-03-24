#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
implements classes for writing RNA secondary structure files
"""

__author__ = "Tomasz Puton"
__credits__ = "Kristian Rother"
__license__ = "GNU GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

from .secstruc import BasePairs

class BpseqWriterError(Exception): pass

class BpseqWriter:
    """Creates a file with RNA secondary structure in bpseq format"""
    
    def __init__(self, bplist, seq, output_path):
        """Initializes BpseqWriter.
        
        bplist: a BasePairs instance
        seq: sequence of an RNA chain
        output_path: path to a new file, where RNA sec struct in bpseq format
            will be written
        """
        if not isinstance(bplist, BasePairs):
            raise BpseqWriterError("%s is not an instance of BasePairs class" \
            % bplist)
            
        self.output_file = open(output_path,'w')
        self.bplist = bplist.directed()
        self.seq = seq
        self.bpseq = ""
        
    def calculateBpseq(self):
        """Calculates RNA sec struct in bpseq format"""
        bpseq = ""
        partners = self.bplist.toPartners(len(self.seq), -1)

        for position, partners, nt in \
            zip(list(range(1, len(self.seq) + 1)), partners, self.seq):

            if partners == None:
                bpseq += str(position) + ' ' + nt + ' ' + '0\n'
            else:
                if len(partners) == 1:
                    bpseq += str(position) + ' ' + nt + ' ' + '%d\n' \
                    % (partners[0] + 1,)
                else:
                    for item in sorted(partners):
                        bpseq += str(position) + ' ' + nt + ' ' + '%d\n' \
                        % (item + 1,)
                        
        self.bpseq = bpseq

    def makeBpseqFile(self):
        """Writes RNA sec struct in bpseq format to a file on a disk"""
        self.output_file.write(self.bpseq)
        self.output_file.close()
        
    
