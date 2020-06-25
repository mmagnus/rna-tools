#!/usr/bin/env python
#
# validators.py
#
# RNA 3D structure prediction.
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

"""Checking of toplevel command parameters."""

from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.Template import Template
from rna_tools.tools.mini_moderna3.moderna.RNAModel import RnaModel
from rna_tools.tools.mini_moderna3.moderna.ModernaFragment import ModernaFragment
from rna_tools.tools.mini_moderna3.moderna.fragment_library.SearchLIR import LirHit, FragmentCandidates
from rna_tools.tools.mini_moderna3.moderna.RNAResidue import RNAResidue
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaAlphabet import alphabet
from rna_tools.tools.mini_moderna3.moderna.util.Errors import ParameterError, AlphabetError
from rna_tools.tools.mini_moderna3.moderna.sequence.RNAAlignment import RNAAlignment, read_alignment
import os, re

def validate_structure(struc):
    """Checks for ModernaStructure objects."""
    if isinstance(struc, ModernaStructure):
        return struc
    raise ParameterError("Bad parameter: '%s' must be a structure object."%str(struc))
    
def validate_template(template):
    """Checks Template object"""
    if isinstance(template, Template):
        return template
    raise ParameterError("Bad parameter: '%s' must be a Template object."%str(template))

def validate_model(model):
    """Checks RnaModel object"""
    if isinstance(model, RnaModel):
        return model
    raise ParameterError("Bad parameter: '%s' must be a RnaModel object."%str(model))

def validate_fragment(frag):
    """Checks RnaModel object"""
    if isinstance(frag, ModernaFragment) or isinstance(frag, LirHit):
        return frag
    raise ParameterError("Bad parameter: '%s' must be a ModernaFragment object."%str(frag))

def validate_alignment(ali):
    """Checks alignment string or object"""
    if isinstance(ali, RNAAlignment):
        return ali
    elif isinstance(ali, str):
        return read_alignment(ali)
    raise ParameterError("Bad parameter: '%s' must be a pairwise Alignment (object or fasta string)."%str(ali))
    
def validate_seq(seq):
    """Checks sequence format"""
    if isinstance(seq, Sequence):
        return seq
    elif type(seq) == str:
        seq = Sequence(seq)
        return seq
    raise ParameterError("Bad parameter: '%s' must be a RNA sequence."%str(seq))

def validate_secstruc(secstruc):
    """Checks vienna string format"""
    if re.search("^[\(\)\.]*$", secstruc):
        return secstruc
    raise ParameterError("Bad parameter: '%s' must be a RNA secondary structure in dot-bracket format."%str(secstruc))

def validate_resnum(resnum):
    """Checks residue number format"""
    resnum = str(resnum)
    if re.search("^\d+\w*$", resnum):
        return resnum
    raise ParameterError("Bad parameter: '%s' must be a residue number."%str(resnum))

def validate_resnum_list(resnum_list, length=None):
    """Checks list of residues"""
    if hasattr(resnum_list, '__iter__'):
        resnum_list = [validate_resnum(resi) for resi in resnum_list]
        if length==None or length==len(resnum_list):
            return resnum_list
    raise ParameterError("Bad parameter: '%s' must be a list of residues."%str(resnum_list))

def validate_frag_candidate_list(fc_list):
    if isinstance(fc_list, FragmentCandidates):
        return fc_list
    raise ParameterError("Bad parameter: '%s' must be a list of fragment candidates."%str(fc_list))
    
def validate_resi(resi):
    """Checks RNAResidue"""
    if isinstance(resi, RNAResidue):
        return resi
    raise ParameterError("Bad parameter: '%s' must be a ModeRNAResidue."%str(resi))
    
def validate_resi_list(resi_list, length=None):
    """Checks list of residues"""
    if hasattr(resi_list, '__iter__'):
        resi_list = [validate_resi(resi) for resi in resi_list]
        if length==None or length==len(resi_list):
            return resi_list
    raise ParameterError("Bad parameter: '%s' must be a list of residues."%str(resi_list))
    
def validate_alphabet(alpha):
    """Checks residue names."""
    if str(alpha) in alphabet:
        return alphabet[str(alpha)].long_abbrev
    else:
        try:
            return alphabet.get_short_original(str(alpha)).long_abbrev
        except AlphabetError:
            pass
    raise ParameterError("Bad parameter: '%s' must be a valid residue name."%str(alpha))

def validate_alphabet_list(alphabet_list, length=None):
    """Checks list of residue names"""
    if hasattr(alphabet_list, '__iter__') or isinstance(alphabet_list, str):
        alphabet_list = [validate_alphabet(alpha) for alpha in alphabet_list]
        if length==None or length==len(alphabet_list):
            return alphabet_list
    raise ParameterError("Bad parameter: '%s' must be a list of residue names."%str(alphabet_list))

    
def validate_filename(filename, exist=False):
    """Checks file names"""
    if exist:
        if not os.path.exists(filename):
            raise ParameterError("Bad parameter: '%s' must be a valid file name."%str(filename))
    else:
        path = os.path.dirname(filename)
        validate_path(path)
    return filename
    
def validate_path(path, exist=True):
    """Checks path names."""
    if exist:
        if path != '' and not os.path.exists(path):
            raise ParameterError("Bad parameter: '%s' must be a valid path name."%str(path))
    return path
        
    
    
    
    
