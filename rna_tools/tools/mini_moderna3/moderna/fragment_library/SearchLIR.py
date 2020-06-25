#!/usr/bin/env python
#
# SearchLIR.py
#
# Searches a local LIR file for RNA fragments that fits.
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

import os
from numpy import array
# from rna_tools.tools.mini_moderna3.moderna.fragment_library.LIR import Lir, LirRecord
from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.ModernaFragment import keep_first_last
from rna_tools.tools.mini_moderna3.moderna.FragmentInsertion import FragmentInserter
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence
from rna_tools.tools.mini_moderna3.moderna.fragment_library.StructureLibrary import library
from rna_tools.tools.mini_moderna3.moderna.util.Errors import LirError, SearchLirError

from rna_tools.tools.mini_moderna3.moderna.Constants import PATH_TO_LIR_STRUCTURES, \
                LIR_DATABASE_PATH, MAX_DIST_STEM, \
                NUMBER_OF_FRAGMENT_CANDIDATES

from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log

struc_cache = {}

DIST_MASK = array([0, 1, 1, 1, 1, 0, 0, 1, 1])

class LirScoringOption(object):
    """
    Class for gathering all Lir scoring values toogether.
    Has an influence on the way in which fragment candidates are choosen. 
    Enables to distinguish between fast scoring functions and advanced scoring functions.
    """
    def __init__(self, scoring_mode = 'fast'): 
        if scoring_mode == 'fast': self.set_fast_scoring()
        elif scoring_mode == 'advanced': self.set_advanced_scoring()


    def set_fast_scoring(self, atom_distance = 10.0, seq_sim=0.0, secstruc=100.0):
        # previously tried: omega 0.0
        """
        Enables to change score for goodness and sequence similarity.
        RMSD, clashes and HBonds are always 0 for fast scoring.
        """    
        self.distance = atom_distance
        self.seq_similarity = seq_sim
        self.rmsd = 0
        self.clashes = 0
        self.secstruc = secstruc
    

    def set_advanced_scoring(self, atom_distance = 0.0000, seq_sim=1.0, rmsd=10.0, clash=2.00, secstruc=100.0):
        # previously: omega 5.0, good 2.0 rms 10.0
        """
        Enables to change all scoring options
        (goodness, seqence similarity, RMSD, clashes, H bonds).
        """   
        self.distance = atom_distance  
        self.seq_similarity = seq_sim
        self.rmsd = rmsd
        self.clashes = clash     
        self.secstruc = secstruc


class LirRecord():
    """empty class @mmagnus"""
    pass

class LirQuery(LirRecord):
    """
    Special LirRecord taht enables searching lfragments.
    """
    def __init__(self, res5, res3, sequence, model_structure, lir_path, secstruc=None):
        LirRecord.__init__(self, fr_length = len(sequence), structure = None, chain = None, preceding_resi = None, following_resi=None, \
            sequence=sequence, sequence_anchor=None, secstruc=secstruc, x=None, y=None, dist_anchor=None, beta=None, gamma=None, omega5=None, omega3=None, \
            P_dist=None, O5p_dist=None, C5p_dist=None, C4p_dist=None, C3p_dist=None, O3p_dist=None, O2p_dist=None, C1p_dist=None, N_dist=None)
        self.anchor5 = res5 # ModernaResiude instance
        self.anchor3 = res3 # ModernaResidues instance
        self.sequence = sequence # ModernaSequence instance/string ???
        
        self.model_instance = model_structure # RNAModel instance
        self.set_record_values() 
        self.lir_path = lir_path
        if secstruc and len(sequence)+2!=len(secstruc):
            raise SearchLirError("Secondary structure '%s' must be 2 resis shorter as the sequence '%s'."%(secstruc, str(sequence)))
#        Whether function get_query_record in LIR should exist?
        

    def set_record_values(self):
        l=Lir(self.anchor5, self.anchor3)
        self.fr_length = len(self.sequence) 
        self.x=l.x
        self.y=l.y
        self.dist_anchor = l.dist_anchor
        self.beta = l.beta
        self.gamma = l.gamma
        self.omega5 =l.omega5
        self.omega3 = l.omega3
        self.distances = array([l.P_dist, l.O5p_dist, l.C5p_dist, l.C4p_dist, l.C3p_dist, l.O3p_dist, l.O2p_dist, l.C1p_dist, l.N_dist])


class LirHit(LirRecord):
    """
    Special LIRRecord that has a goodness, rmsd and sequence identity score.
    """

    def __init__(self, query, fr_length, structure, chain, preceding_resi, following_resi, \
        sequence, sequence_anchor, secstruc, x, y, dist_anchor, beta, gamma, omega5, omega3,  \
        P_dist, O5p_dist, C5p_dist, C4p_dist, C3p_dist, O3p_dist, O2p_dist, C1p_dist, N_dist):
        LirRecord.__init__(self, fr_length, structure, chain, preceding_resi,following_resi, \
            sequence, sequence_anchor, secstruc, x, y, dist_anchor, beta, gamma, omega5, omega3, 
            P_dist, O5p_dist, C5p_dist, C4p_dist, C3p_dist, O3p_dist, O2p_dist, C1p_dist, N_dist)
        
        self.d_dist_score = 0
        self.rmsd = 0
        self.seq_similarity = 0
        self.clash = 0 #score depends on the number of clashing residues
        self.score = 0
        self.d_secstruc = 0
        self.fragment_instance = None
        self.query = query
    
    def __cmp__(self, hit):
        """Allow lists of LirHits to be sorted."""
        return cmp(self.score, hit.score)


    def __str__(self):
        """Returns one-line summary."""
        hit_str="%16s [%s]\t%4s\t%s\t%s\t%3i\t%7.4f\t%7.4f\t%7.4f\t%7.4f"\
                % (self.structure,self.chain,self.preceding_residue, self.sequence, self.secstruc, self.fr_length,self.d_dist_score, self.rmsd, self.clash, self.score)
        return hit_str

    def __repr__(self):
        return str(self.score)

    def calculate_dist_score(self):
        """calculates square differences of all distances."""
        self.d_dist_score = sum( (self.distances*DIST_MASK - self.query.distances*DIST_MASK)**2)


    def calculate_seq_similarity(self):
        """
        The higher score the less similar sequences are.
        Two identical sequences will get score 0.
        """
        #KR: used local variables to save time
        hit_seq = self.sequence
        query_seq = self.query.sequence
        assert len(hit_seq) == len(query_seq)
        # calculate similarity
        similarity_score = 0.0
        for hit_ae, query_ae in zip(hit_seq, query_seq):
            if hit_ae == query_ae: continue
            hit_orig = hit_ae.original_base
            query_orig = query_ae.original_base
            if hit_orig == query_orig: similarity_score += 0.1
            elif hit_orig in 'AG' and query_orig in 'AG': similarity_score += 0.5
            elif hit_orig in 'CU' and query_orig in 'CU': similarity_score += 0.5
            else: similarity_score +=1
                    
        if len(self.sequence): self.seq_similarity = similarity_score/len(self.sequence)
        else: self.seq_similarity = 0

    def match_secstruc(self):
        """Gives a BIG penalty if secstruc does not match."""
        if self.query.secstruc:
            self.d_secstruc = 100.0
            if self.fragment_instance:
                if self.query.secstruc == self.fragment_instance.struc.get_secstruc():
                    self.d_secstruc = 0.0
            elif self.query.secstruc == self.secstruc:
                self.d_secstruc = 0.0
            
    def get_fragment(self):
        """
        Allows to get small structure (as ModernaFragment) with fragment candidate fragment (fragment + anchor residues).
        Returns ModernaFragment. Fragment is also stored in fragment_instance attribute.

        Argument:
        - anchor residue from 5' end (from structure to which fragment will be added)
        - anchor residue from 3' end (from structure to which fragment will be added)
        - sequence
        """
        if not self.fragment_instance:
            lir_path = self.query.lir_path
            frag = library.get_fragment_part(\
                    lir_path+self.structure, self.chain, \
                    self.preceding_residue, self.following_residue, \
                    self.query.anchor5, self.query.anchor3, \
                    self.query.sequence, keep=keep_first_last, \
                    seq=self.sequence_anchor \
                    )
            finsert = FragmentInserter()
            finsert.prepare_fragment(frag, self.query.model_instance)
            self.fragment_instance = frag
        return self.fragment_instance

    def calculate_rmsd(self, res5=None, res3=None, seq=None):
        """
        Fragnet is superimposed acording coordinates of given anchor residues (taken from query) and RMSD is calculated.
        Returns RMSD. RMSD is also stored in rmsd attrbute of fragment hit.
        """
        self.get_fragment() 
        self.rmsd = self.fragment_instance.rmsd
        return self.rmsd

    def find_clash(self):
        """
        Checks wheader a fragment candidate clashes with given list of residues or not.
        Returns list with tuples of clashing residues.
        
        Arguments:
        - list with residues from structure for which fragment will be added
        - anchor residue from 5' end (from structure to which fragment will be added)
        - anchor residue from 3' end (from structure to which fragment will be added)
        - fragment sequence (from structure to which fragment will be added)
        """
        self.get_fragment()
        model_residues = self.query.model_instance.find_residues_not_in_range(self.query.anchor5.identifier,  self.query.anchor3.identifier)
        clashes = self.fragment_instance.has_clashes(model_residues)
        self.clash = clashes
        return clashes
    
    def score_clash(self):
        """
        Finds clashes and gives score acording to number of clashes.
        If there is no clash the score is 0.
        """
        clashes = self.find_clash()
        self.clash = len(clashes)
        return self.clash
    
    def check_backbone_continouity(self): 
        pass
    
    def calculate_score(self,  scoring):
        """
        Calculates hit score acording to the provided scoring_option.
        Scoring_option instance is responsible for distingushing fast scoring from advance scoring
        The lower score the better candidate.
        
        Atributes:
        - LirScoringOption instance
        """
        if scoring.distance: self.calculate_dist_score()
        if scoring.seq_similarity:  self.calculate_seq_similarity()
        if scoring.rmsd: self.calculate_rmsd()
        if scoring.clashes:  self.score_clash()
        if scoring.secstruc: self.match_secstruc()
        
        self.score = scoring.distance * self.d_dist_score + \
                     scoring.seq_similarity * self.seq_similarity + \
                     scoring.rmsd * self.rmsd + \
                     scoring.clashes * self.clash + \
                     scoring.secstruc * self.d_secstruc
    

    def write_hit_structure(self,  file_name='LirHit.pdb', with_anchor_residues=False, with_model=False,  write_anchors_to_file=True):
        """
        Writes structure of hit to a pdb file.
        Optioonaly it can write hit with anchor residues or with whole model.
        
        Arguments:
        - with_anchor_residues - True/False (by default False)
        - with_model - True/False (by default False)
        """
        
        if not self.fragment_instance: self.get_fragment()
        if not self.rmsd: self.fragment_instance.superimpose()   
        
        if write_anchors_to_file: self.write_anchors_to_separate_file(file_name)        
        if with_anchor_residues and not with_model: 
            resis = list(self.fragment_instance.struc)
        elif not with_anchor_residues and not with_model: 
            #with model or without model and without anchor residues
            resis = self.fragment_instance.nonanchor_residues
        if with_model:
            # get a variable with the class because it cant be imported; looks ugly but works.
            rm_class = self.query.model_instance.__class__
            m = rm_class(None, None, self.query.model_instance.chain_name, 'residues', self.query.model_instance)
            m.insert_fragment(self.fragment_instance) 
            m.fix_backbone()
            resis = m.moderna_residues.values() # MM: or resi = m.moderna_residues
                    
        m = ModernaStructure('residues',  resis,  'A')
        m.sort_residues()
        m.write_pdb_file(file_name)
    
    
class FragmentCandidates(object): 
    """
    Takes care about fragment candidates.
    Prepares initial candidates list and different scorings on it.

    Arguments:
    - query (LirQuery instance)
    - path to file with LIR data   
    """
    lir_cache = ["", []] # filename, records
    
    def __init__(self, query, lir_db_filename=LIR_DATABASE_PATH):
        self.query=query
        self.lir_database_filename = lir_db_filename 
        self.accepted_fragments=[]
        self.index = 0


    def __getitem__(self, args):
        return self.accepted_fragments[args]

    
    def __len__(self):
        return len(self.accepted_fragments)
        
    
    def __iter__(self):
        return self.accepted_fragments.__iter__()
 
    def __str__(self):
        return "Present number of accepted fragment candidates: %s" %str(len(self.accepted_fragments))

    def parse_lir_database(self, separator='\t'): # WAS: ' - '
        """
        Reads fragments from the LIR database file and generates 
        LirRecords from the columns of the file.

        To lern more see also LirRecord documentation.
        """
        # check if a different file is read
        if self.lir_cache[0] != self.lir_database_filename or self.lir_cache[1]==[]:
            # empty old lir
            self.lir_cache[1] = []
            # read new lir
            self.lir_cache[0] = self.lir_database_filename
            
            for line in open(self.lir_database_filename):
                if line[0] == 'l': continue #.startswith('fr length'): continue
                line = line.strip().split(separator)
                if len(line) == 17:
                    self.lir_cache[1].append(line)
                    
    def create_hit_from_line(self, line):
        """Creates a LirHit object."""
        lir_hit = LirHit(
            query = self.query, 
            fr_length=int(line[0]), 
            structure=line[1], 
            chain=line[2],
            preceding_resi=line[3],
            following_resi=line[4], 
            sequence = Sequence(line[5]),
            sequence_anchor=Sequence(line[6]), 
            secstruc = line[7], 
            x=0.0,  #float(line[7]),
            y=0.0 ,#float(line[8]),
            dist_anchor =0.0,  #float(line[9]),
            beta =0.0,  #float(line[10]),
            gamma=0.0,  #float(line[11]), 
            omega5=0.0,#float(line[12]), 
            omega3=0.0,#float(line[13])
            P_dist = float(line[8]), 
            O5p_dist = float(line[9]), 
            C5p_dist = float(line[10]), 
            C4p_dist = float(line[11]), 
            C3p_dist = float(line[12]), 
            O3p_dist = float(line[13]), 
            O2p_dist = float(line[14]), 
            C1p_dist = float(line[15]), 
            N_dist = float(line[16])
        )
        return lir_hit


    def create_initial_fragment_set(self,  max_dist_anchor=MAX_DIST_STEM,  separator='\t'):
        """
        Collects all records from the LIR database that fill the length and anchor distance conditions.
        """
        self.parse_lir_database(separator)
        #Creates a list of LirHit objects for all entries in the pre-parsed LIR DB.
        length = self.query.fr_length # SPEEDUP
        self.accepted_fragments = [self.create_hit_from_line(line) for line in self.lir_cache[1] if int(line[0])==length]


    def make_fast_scoring(self,  scoring_option,  number_of_candidates):
        """
        Prepares scoring for all candidates based on values that are already present in LIR database.
        Takes into account goodness ans sequence similarity.
        Returns small number of candidates (indicated by number of candidates) 
        which can udergo more advance (longer) scoring.
        """
        #for hit in self.accepted_fragments: hit.calculate_score(scoring_option)
        [hit.calculate_score(scoring_option) for hit in self.accepted_fragments]
        self.accepted_fragments.sort()
        self.accepted_fragments = self.accepted_fragments[:number_of_candidates]
  
  
    def make_advanced_scoring(self, scoring_option):
        """
        Prepares scoring acorgding to RMSD, clashes and Hbonds.
        Should be called when number of self.accepted_fragments is small 
        so as to avoid long calculations.
        """
        #for hit in self.accepted_fragments: hit.calculate_score(scoring_option)
        [hit.calculate_score(scoring_option) for hit in self.accepted_fragments]
        self.accepted_fragments.sort()
        if self.accepted_fragments != []:
            if not self.accepted_fragments[0].fragment_instance: 
                for hit in self.accepted_fragments: hit.get_fragment()
        #TODO: could be merged with method above.


    def write_fragment_candidates(self, directory_name='LIR_candidates',  with_anchor_residues=False,  with_model=False, write_anchors_to_file=False,  log=True):
        """
        Writes pdb files with all fragment candidates.
        Candidates can be easly checked manualy by calling 'pymol *' in candidates directory
        
        Arguments:
        - directory name 
        """
        # MM: What when there is no hit? Exception or nothing?
        if not os.path.isdir(directory_name): 
            os.mkdir(directory_name)
        if directory_name[-1] != os.sep: directory_name += os.sep
        if log: f=open(directory_name+'candidates.log', 'w')
        for x, hit in enumerate(self.accepted_fragments):
            if log: f.write(str(hit)+'\n'+50*'_'+'\n')
            path=directory_name  + 'candidate' + str(x) + '.pdb'
            hit.write_hit_structure(path, with_anchor_residues,  with_model,  write_anchors_to_file)



class FragmentFinder(object):
    """
    Looks for fragment that fits the best."""
    def __init__(self,  anchor5, anchor3, sequence, model_structure,  candidates_number=NUMBER_OF_FRAGMENT_CANDIDATES, lir_path=PATH_TO_LIR_STRUCTURES, secstruc=None):
        self.anchor5 = anchor5 # ModernaResidue instance, anchor residue from the model from 5' (fragment start) 
        self.anchor3 = anchor3 # ModernaResidue instance, anchor residue from the model from 3' (fragment end)
        self.sequence = sequence # Sequence instance. The sequence that fragment should have.
        self.model_structure = model_structure #RNAModel instance
        self.candidates_number = candidates_number
        self.secstruc = self._get_secstruc(secstruc)
        self.scoring = LirScoringOption()
        self.lir_path=lir_path
        self.query = self.get_query()
        #TODO: could check whether parameters have the right type (?)
        #MM: but type(anchor5) gives  only 'instance'
    
    def _get_secstruc(self, secstruc):
        """Adds pairing of anchors to secstruc query"""
        if secstruc == None:
            return None
        basepair = self.anchor5.get_bp(self.anchor3)
        if basepair and (basepair.canonical or basepair.wobble):
                return '('+secstruc+')'
        else:
            return '.'+secstruc+'.'
            
    def get_query(self):
        """Creates LirQuery obiect taht can be used for searching LIRdb."""
        return LirQuery(self.anchor5, self.anchor3, self.sequence, self.model_structure, self.lir_path, self.secstruc)

    def log_candidates(self, candidates):
        log.write_message('\nFragment candidates:\n') 
        log.write_message("""
Fragment candidates:
structure [chain]\t5'-resi\tsequence\tsecstruc\tlength\tdist\tRMSD\tclashes\tscore""")
        for l in candidates: log.write_message(str(l))

    def find_fragment_candidates(self):
        """Returns a FragmentCandidates object."""
        candidates = FragmentCandidates(self.query)
        candidates.create_initial_fragment_set()        
        self.scoring.set_fast_scoring()
        candidates.make_fast_scoring(self.scoring, self.candidates_number)
        self.scoring.set_advanced_scoring()
        candidates.make_advanced_scoring(self.scoring)
        self.log_candidates(candidates)        
        if len(candidates)>0 and candidates[0].score > 10000:
            log.write_message("\nNo fragment candidate with aproppriate secondary structure was found.\nThe next best fragment is inserted.\n")
        return candidates
                
    def find_fragment(self):
        """
        Search for fragment candidates and returns fragment that belongs to the best candidate.
        """
        candidates = self.find_fragment_candidates() # FragmentCandidates instance
        if len(candidates) == 0: raise LirError('No fragment candidates found')
        return candidates[0].fragment_instance
