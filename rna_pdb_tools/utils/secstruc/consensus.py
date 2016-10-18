#!/usr/bin/env python
#!-*-coding: utf-8-*-

"""
Tools for calculating consensus RNA secondary structure.
"""

__author__ = "Tomasz Puton"
__license__ = "GNU GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

import numpy as np
from secstruc import ViennaStructure, BasePairs

class ConsensusSecstrucError(Exception): pass

__all__ = ['calc_consensus_secstruc_for_bps',
           'calc_consensus_secstruc_for_dbn']


def get_length(secstrucs):
    """Checks that all secondary structures in the set have the same length."""
    n = len(secstrucs[0])
    for s in secstrucs:
        if len(s) != n:
            raise ConsensusSecstrucError("Length mismatch in the secondary structures: %s"""%str(secstrucs))
    return n


def make_matrix_with_indexes_as_cells(dim1, dim2):
    """Generates numpy array with tuples as strings representing indexes of
        paired residues.
        
        dim1: int, first dimention of a matrix
        dim2: int, second dimention of a matrix
    """
    empty = np.empty(dim1*dim2, dtype='|S15').reshape(dim1, dim2)
    for i in range(dim1):
        for j in range(dim2):
            empty[i,j] = "(%i,%i)" % (i,j)
    return empty


def turn_bps_into_interaction_matrix(bps, length, value=1, offset=-1):
    """Converts RNA secondary structure in BasePairs object into contact map
        (numpy array), where 0 means no pair, 1 pair.
        
        Important! 
        
        bps: BasePairs object
        length: int, length of an RNA sequence
        value: 1 or float, value to be placed in the appropriate matrix cell
            in order to denote a pair
        offset: int, compensates the problem of numbering residues.
            Should be set to -1, if numbering of bases in BasePairs
            starts from 0, not from 1!!! This is the case when input is
            BasePairs object generated from the output of RNAView, i.e. by
            SecondaryStructureCalculator.get_base_pairs !!!
            
        It returns a matrix of residue pairs generated using numpy.
    """
    matrix = np.zeros(length * length).reshape(length, length)
    for row, col in bps:
        # the matrix does not need to be symmetrical
        matrix[row+offset, col+offset] = value
    return matrix


def turn_dbn_into_interaction_matrix(dbn, value=1, offset=0):
    """Converts RNA secondary structure in Vienna format (dot-bracket notation)
        into contact map (numpy array), where 0 means no pair, 1 pair.
        
        dbn: str, RNA secondary structure in Vienna format e.g. ...((...))..
        value: 1 or float, value to be placed in the appropriate matrix cell
            in order to denote a pair
            
        It returns a matrix of residue pairs generated using numpy.
    """
    len_dbn = len(dbn)
    pairs = ViennaStructure(dbn).toPairs(offset=offset)
    return turn_bps_into_interaction_matrix(pairs, len_dbn, value=value,
                                            offset=offset)


def get_consensus_base_pairs(secstrucs, length=None, threshold=1.0):
    """Returns consensus BasePairs object for either a list of strings with
        RNA secondary structure in Vienna format or list of BasePairs objects
        ('secstruc' argument).
        
        secstrucs: a list of strings with RNA secondary structure in Vienna
            format or list of BasePairs objects
        length: optional int; only used if secstrucs contains BasePairs objects
        threshold: float, threshold for calculating consensus 
    """
    if len(secstrucs) == 0:
        raise ConsensusSecstrucError("No data for calculating consensus!")
        
    input_type = set(type(x) for x in secstrucs)
    input_type = tuple(input_type)
    if len(input_type) != 1:
        error_msg = "Expected only one type, got %s" % str(input_type)
        raise ConsensusSecstrucError(error_msg)
    
    input_type = input_type[0]
    if input_type is str:
        input_lengths = [len(x) for x in secstrucs]
        input_lengths = tuple(set(input_lengths))
        
        if len(input_lengths) != 1:
            error_msg = "Input Viennas have different length, "+\
                        "got %s" % str(input_lengths)
            raise ConsensusSecstrucError(error_msg)
            
        length = input_lengths[0]
        
        if length is not None and input_lengths[0] != length:
            error_msg = "Length (%d) was specified manually, " % length +\
                        "but it differs from length of "+\
                        "input Viennas (%d)." % input_lengths[0]
            raise ConsensusSecstrucError(error_msg)
            
        matrix_generator = turn_dbn_into_interaction_matrix
        kwargs = {'offset' : 0}
            
    elif input_type is BasePairs:
        if length is None:
            raise ConsensusSecstrucError("length kwarg is compulsory!")
            
        matrix_generator = turn_bps_into_interaction_matrix
        kwargs = {'length' : length, 'offset' : -1}
        
    else:
        error_msg = "secstrucs should contain strings or BasePairs objects"
        raise ConsensusSecstrucError(error_msg)

    if not 0.0 < threshold <= 1.0:
        err_msg = "Invalid threshold value %4.2f. Expected from range (0, 1>" \
                  % threshold
        raise ConsensusSecstrucError(err_msg)

    matrix = matrix_generator(secstrucs[0], **kwargs)
    for secstruc in secstrucs[1:]:
        matrix += matrix_generator(secstruc, **kwargs)
        
    matrix = matrix / len(secstrucs)
    consensus_pairs = matrix >= threshold
        
    indexes = make_matrix_with_indexes_as_cells(length, length)
    # tu trzeba uwzględnić kwargs['offset']
    
    bps = BasePairs(eval(x) for x in indexes[consensus_pairs])
    return BasePairs((x[0]-kwargs['offset'], x[1]-kwargs['offset']) for x in
        bps)


def calc_consensus_secstruc_for_dbn(secstrucs, threshold=1.0):
    """
    Consensus Secondary Structure Calculator
    Takes a list of secondary structure strings in the
    dot-bracket format, and calculates a consensus structure from them.
    
        secstrucs: list of strings, each string is RNA secondary structure in
            Vienna format
        threshold: float, threshold for calculating consensus
    
    Returns consensus Vienna.
    """
    bps = get_consensus_base_pairs(secstrucs, threshold=threshold)
    return bps.toVienna(len(secstrucs[0]), offset=0)
    
    
def calc_consensus_secstruc_for_bps(bps, length=None, threshold=1.0):
    """
    Consensus Secondary Structure Calculator
    Takes a list of BasePairs objects and calculates a consensus structure
    from them.
    
        bps: list of BasePairs objects
        length: length of an RNA for which BasePairs come from
        threshold: float, threshold for calculating consensus
    
    Returns consensus BasePairs object.
    """
    if length is None:
        raise ValueError("length cannot be None")
    return get_consensus_base_pairs(bps, length, threshold)


def _most_counted_char(count, min_count, i):
    """find the most counted element"""
    rev_items = [(count[k],k) for k in count]
    if len(rev_items) > 0:
        rev_items.sort()
        max_item = rev_items[-1]
        if max_item[0] >= min_count:
            bp = max_item[1]
            if bp > i:
                return '('
            else:
                return ')'
    return '.'
