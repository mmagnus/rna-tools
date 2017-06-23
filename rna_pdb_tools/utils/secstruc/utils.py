#!/usr/bin/env python
#-*-coding: utf-8-*-

"""
Various helper functions for processing RNA secondary structure.
"""

__author__ = "Tomasz Puton"
__credits__ = "Kristian Rother"
__license__ = "GNU GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Tomasz Puton"
__email__ = "t.puton@amu.edu.pl"
__status__ = "Production"

import os, tempfile, re, copy, types
from math import sqrt

# @magnus hack!!!!
#from wrappers.alignment_local import RCoffee, CMalign
#from wrappers.paralign_rfamseq import ParalignRfamseq, \
#ParalignTableParserRfamseqError
from .secstruc import BasePairs

from cogent import LoadSeqs
from cogent.struct.pairs_util import check_structures, extract_seqs, \
selectivity_formula, sensitivity_formula
from cogent.core.alignment import Alignment, SequenceCollection

from cogent.struct.pairs_util import ACCEPTED as ACCEPTED_STD
# ACCEPTED_EXT contains possible base pairs from 'extended secondary structure'.
# Basically, all possible base pairs are included in this list.
ACCEPTED_EXT = dict.fromkeys(list(map(tuple, ["AA", "AC", "AG", "AU", "CA", "CC",
                                         "CG", "CU", "GA", "GC", "GG", "GU",
                                         "UA", "UC", "UG", "UU"])))

def mcc_formula(counts):
    """Return correlation coefficient
    
    counts: dict of counts, containing at least TP, TN, FN, FP, and FP_COMP
    """
    tp = counts['TP']
    tn = counts['TN']
    fp = counts['FP']
    fn = counts['FN']
    comp = counts['FP_COMP']
    
    mcc_quotient = (tp+fp-comp)*(tp+fn)*(tn+fp-comp)*(tn+fn)
    if mcc_quotient > 0:
        mcc = (tp*tn-(fp-comp)*fn)/sqrt(mcc_quotient)
    else:
        #raise ValueError("mcc_quotient <= 0: %.2f"%(mcc_quotient))
        mcc = 0 # @hack @magnus
    return mcc

def get_all_pairs(sequences, min_dist=1, accepted=ACCEPTED_STD):
    """Return number of possible base pairs in the sequece
    
    sequences: list of Sequence objects or strings
    min_dist: integer, minimum distance between two members of a 
        base pair. Default is 4 (i.e. minimum of 3 unpaired bases in a loop)

    The number of pairs is defined as all possible GC, AU, and GU pairs, 
        respecting the minimum distance between the two members of a
        base pair.
    This method returns the average number of possible base pairs over
        all provided sequences.
    """
    # TAKEN FROM cogent.struct.pairs_utils.py but refactored:
    # added third argument accepted, to be passed to funcion.
    if min_dist < 1:
        raise ValueError("Minimum distance should be >= 1")
    if not sequences:
        return 0.0
    tn_counts = []
    for seq in sequences:
        seq_str = str(seq).upper()
        seq_count = 0
        #print 'xrange', range(len(seq)-min_dist)
        for x in range(len(seq)-min_dist):
            for y in range(x+min_dist,len(seq)):
                if (seq_str[x], seq_str[y]) in accepted:
                    #print x,y, seq_str[x], seq_str[y], 'Y'
                    seq_count += 1
                else:
                    pass
                    #print x,y, seq_str[x], seq_str[y], 'N'
        tn_counts.append(seq_count)
    return sum(tn_counts)/len(tn_counts)
    
def get_counts(ref, predicted, split_fp=True, sequences=None, min_dist=1,
               accepted=ACCEPTED_STD):
    """Return TP, TN, FPcont, FPconf FPcomp, FN counts
    
        ref: BasePairs object
        pred: BasePairs object
        split_fp: bool, set flags whether to split false positivies into
            different categories or not
        sequences: list of sequences of RNA (strings)
        min_dist: minimum distance between two bases which make up a pair
        accepted: dict, pairs of letters denoting allowed base pairs
    """
    assert isinstance(ref, BasePairs), 'ref is not BasePairs instance'
    assert isinstance(predicted, BasePairs), 'predicted is not BasePairs instance'
    assert isinstance(split_fp, bool), 'split_fp is not bool'
    assert type(sequences) in [list, tuple, type(None)],\
        'sequences is not tuple or list'
    assert isinstance(min_dist, int), 'min_dist is not integer'
    assert isinstance(accepted, dict), 'accepted is not dict'
    
    from cogent.struct.pairs_util import get_counts
    counts = get_counts(ref=ref, predicted=predicted, split_fp=split_fp,
                        sequences=None, min_dist=min_dist)
    
    if sequences:
        num_possible_pairs = get_all_pairs(sequences, min_dist,
                                           accepted=accepted)
        counts['TN'] = num_possible_pairs - counts['TP'] -\
            counts['FP_INCONS'] - counts['FP_CONTRA']
        
    return counts
    
def sensitivity(ref, predicted, accepted=ACCEPTED_STD):
    """Return sensitivity of the predicted structure

    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure
    accepted: dict, pairs of letters denoting allowed base pairs

    Formula: sensitivity = tp/(tp + fn)
    tp = True positives
    fn = False negatives
    """
    check_structures(ref, predicted)
    if not ref and not predicted:
        return 1.0
    elif not predicted:
        return 0.0

    counts = get_counts(ref, predicted, accepted=accepted)
    return sensitivity_formula(counts)

def selectivity(ref, predicted, accepted=ACCEPTED_STD):
    """Return selectivity of the predicted structure
    
    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure
    accepted: dict, pairs of letters denoting allowed base pairs

    Formula: selectivity = tp/(tp+fp-fp_comp)
    tp = True positives
    fp = False positives
    fp_comp = compatible fp pairs
    """
    check_structures(ref, predicted)
    if not ref and not predicted:
        return 1.0
    elif not predicted:
        return 0.0

    counts = get_counts(ref, predicted, split_fp=True, accepted=accepted)
    return selectivity_formula(counts)

def mcc(ref, predicted, seqs, min_dist=1, accepted=ACCEPTED_STD):
    """Return the Matthews correlation coefficient

    ref: Pairs object -> reference structure (true structure)
    predicted: Pairs object -> predicted structure
    seqs: list of sequences, necessary to compute the number of true
        negatives. See documentation of extract_seqs function for 
        accepted formats.
    min_dist: minimum distance required between two members of a base pair.
        Needed to calculate the number of true negatives.
    accepted: dict, pairs of letters denoting allowed base pairs
    """
    check_structures(ref, predicted)
    if not ref and not predicted:
        return 1.0
    elif not predicted:
        return 0.0
    elif not seqs:
        raise ValueError('No sequence provided!')

    sequences = extract_seqs(seqs)
    counts = get_counts(ref, predicted, sequences=sequences, split_fp=True,\
        min_dist=min_dist, accepted=accepted)
    
    result = mcc_formula(counts)
    if result < -1 or result > 1:
        raise ValueError("mcc not in range <-1, 1>: %.2f"%(result))
    return result 


class CannotMakeSequenceAlignmentError(Exception): pass

MAX_NUMBER_SEQS = 3 # default for building RNA sequence alignment for automated
                    # testing of RNA secondary structure prediction programs
                    
E_VALUE_CUTOFF = '1e-5' # e-value for filtering out random RNA sequences

def make_sequence_alignment(collection, seq_db_path, cm=None,
                            max_nr_seqs=MAX_NUMBER_SEQS):
    """Based on input sequences creates a sequence alignment. If one sequence is
    provided runs Paralign on a non-redundant copy of RfamSeq10 database
    to get other sequences and subsequently aligns them with RCoffee.
    If collection consists of more than one sequence, just
    runs RCoffee and produces an alignment.
    
    collection: an instance of PyCogent's SequenceCollection consisting of
        either one or more sequences or a PyCogent's Alignment instance
    seq_db_path: str, a path to formatted sequence database in FASTA format
    cm: str, a covariance model of RNA
    """
    if isinstance(collection, Alignment):
        return collection
    elif isinstance(collection, SequenceCollection):
        if len(collection.Seqs) == 1:
            # We have to run Paralign, to get up to 20 seqs from RfamSeq10 db.
            fasta = collection.Seqs[0].toFasta() + '\n'
            
            paralign = ParalignRfamseq(str(collection.Seqs[0]), seq_db_path)
            try:
                results = paralign.run()
            except ParalignTableParserRfamseqError as error:
                raise CannotMakeSequenceAlignmentError(str(error))
            counter = 0
            
            for hit in results['hits']:
                if  counter == max_nr_seqs:
                    break
                if hit.expect <= float(E_VALUE_CUTOFF):
                    try:
                        fasta += hit.fasta.replace('T','U')
                    except ValueError as error:
                        msg = \
                        "Problem getting sequence for '%s': ValueError: %s" % \
                        (str(error), hit.seq_id)
                        raise CannotMakeSequenceAlignmentError(msg)
                    counter += 1
        else:
            # Great, no sequences need to be searched for!
            fasta = str(collection)#.toFasta()
        
        if len(LoadSeqs(data=fasta, aligned=False).Seqs) <= 1:
            msg = \
            "Paralign failed to find homologous sequences using cutoff %s" % \
            (E_VALUE_CUTOFF,)
            raise CannotMakeSequenceAlignmentError(msg)
        
        if cm:
            fd, path = tempfile.mkstemp()
            open(path, 'w').write(cm)
            cmalign = CMalign(fasta)
            cmalign.add_cmfile(path)
            results = cmalign.run()
            os.close(fd)
            os.unlink(path)
        else:
            rcoffee = RCoffee(fasta)
            results = rcoffee.run()
        alignment = LoadSeqs(data=results['alignment'])
        return alignment
    else:
        msg = "'collection' is neither Alignment nor SequenceCollection!"
        raise CannotMakeSequenceAlignmentError(msg)
        
class ExtractQuerySeqPredError(Exception): pass

def extract_query_seq_pred(result, sequence, name):
    """Extracts prediction for a query sequence from a
        multiple vienna record (key 'vienna_sec_structs').
    """
    new_result = copy.deepcopy(result)
    if new_result.get('best_predicted_ss') is not None:
        if new_result.get('all_energies') is None:
            new_result['best_predicted_energy'] = None
        
        return new_result
        
    else:
        seq_name_pattern = re.compile(r'>(.+)+\n')
        seq_pattern      = re.compile(r'\n(\w+)\n')
        dbn_pattern      = re.compile(r'\n([\.\(\)\[\]\}\{]+)')
        
        seq_names = seq_name_pattern.findall(result['vienna_sec_structs'])
        seqs      = seq_pattern.findall(result['vienna_sec_structs'])
        dbns      = dbn_pattern.findall(result['vienna_sec_structs'])
        
        try:
            seq_index  = seqs.index(sequence)
        except ValueError as error:
            raise ExtractQuerySeqPredError(str(error))
        else:
            if seq_names[seq_index].find(name) == -1:
                msg = "Alignment seq '%s' doesn't match its expected name '%s'." % \
                      (sequence, name)
                raise ExtractQuerySeqPredError(msg)
            # Not an elegant way of copying dict, but we have to do it this way
            # because copy.deepcopy fails with copying cogent's Sequence object
            #_sequence_ = result['sequence']
            #del result['sequence']
            
            #new_result['sequence'] = _sequence_
            #result['sequence'] = _sequence_
            new_result['best_predicted_ss'] = dbns[seq_index]
            if len(new_result['all_energies']) == 0:
                new_result['best_predicted_energy'] = None
            else:
                new_result['best_predicted_energy'] = \
                    new_result['all_energies'][seq_index]
            del new_result['vienna_sec_structs']
            del new_result['all_energies']
            return new_result
    

def get_rfam_key(rfam_id, rfam_file='rfam.txt'):
    """
    Parses lines of the pattern:
1	1	RF00001	5S_rRNA	5S ribosomal RNA...
    """
    for line in open(rfam_file):
        if re.search('^\d+\t\d+\tRF\d\d\d\d\d\t', line):
            tokens = line.split('\t')
            if tokens[2] == rfam_id:
                return tokens[0], tokens[3], tokens[4]
                
def get_pdb_for_rfam_key(rfam_key, pdb_rfam_key='pdb_rfam_reg.txt'):
    """Returns (PDB_ID, Chain) tuples for the given DB ID."""
    result = []
    for line in open(pdb_rfam_key):
        tokens = line.split('\t')
        if tokens[1] == rfam_key:
            result.append((tokens[2], tokens[3]))
    return result
        

def get_pdb_for_rfam(rfam_id):
    """Generates (PDB_ID, chain) tuples for a given RFAM ID."""
    rfam_key = get_rfam_key(rfam_id)[0]
    pdb_ids = get_pdb_for_rfam_key(rfam_key)
    return pdb_ids


def get_bps_for_aligned_seq(aligned_seq, bps, first_index=0):
    """Extracts from an aligned sequence base pairs which map to that sequence.
        e.g.
        ACUAGCUG-----ACUGA
        BasePairs((2, 17))
        
        will return BasePairs mapped to the unaligned sequence:
        BasePairs((2, 12)
        
        aligned_seq: str, aligned sequence e.g. 'ACUGACUAGC---ACGUACGU'
        bps: BasePairs instance
        first_index: int, sets the starting number i.e. usually 0 or 1
    """
    assert isinstance(aligned_seq, str)
    assert isinstance(bps, BasePairs)
    
    partners = bps.toPartners(len(aligned_seq) + first_index)
    new_base_pairs = set()
    
    for i, nt, pos in zip(list(range(first_index, len(aligned_seq) + first_index)),
                          aligned_seq, partners[first_index:]):
        if nt == '-' or pos is None:
            continue
        
        for partner in pos:
            partner_updated_pos = partner - first_index
            if i in partners[partner] and \
            aligned_seq[partner_updated_pos] != '-':
                
                nr_gaps_before_i = aligned_seq[:i].count('-')
                nr_gaps_from_i_to_partner = \
                    aligned_seq[i:partner_updated_pos].count('-')
                
                if i > partner:
                    continue
                
                total_gaps_to_partner = \
                    nr_gaps_before_i + nr_gaps_from_i_to_partner
                
                new_base_pairs.add(tuple(
                    sorted((i - nr_gaps_before_i,
                            partner - total_gaps_to_partner))))
                
    return BasePairs(new_base_pairs).directed()


