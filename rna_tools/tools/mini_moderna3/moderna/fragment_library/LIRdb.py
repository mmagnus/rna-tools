#!/usr/bin/env python
#
# LIRdb.py
#
# Script for preparing file with information about all possible fragments from given structures.
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


import sys, re, os, os.path

from rna_tools.tools.mini_moderna3.moderna.ModernaStructure import ModernaStructure
from rna_tools.tools.mini_moderna3.moderna.sequence.ModernaSequence import Sequence

from LIR import Lir, LirRecord
from rna_tools.tools.mini_moderna3.moderna.Constants import PATH_TO_LIR_STRUCTURES, PATH_TO_LIR_LIST, ANY_RESIDUE,  \
    UNKNOWN_RESIDUE_SHORT, MAX_FRAGMENT_LENGTH
from rna_tools.tools.mini_moderna3.moderna.util.LogFile import log

# constants to create LIR test DB.
#PATH_TO_LIR_STRUCTURES = '../test/test_data/lir_test_files/'
#PATH_TO_LIR_LIST = '../test/test_data/lir_chains.txt'

class ResidueList:
    """
    Manages a list of residues from which linkers can be generated.
    
    If the chain is discontinuous or contains unknown residues it is cut.
    Eg. ACUAXGAG_UUC ---> [[<A><C><U><A>],[<G><A><G>],[<U><U><C>]]
    """
    #KR: introduced this class to have secstruc connected to the residue list.
    def __init__(self, struc):
        """Creates a residue list from an RNA structure."""
        # remove unknown resis
        for res in struc:
            if res.long_abbrev == UNKNOWN_RESIDUE_SHORT:
                struc.remove_residue(res.identifier)
        # identify separate portions of the chain
        seq = struc.get_sequence().seq_with_modifications
        secstruc = struc.get_secstruc()
        resi_list = [res for res in struc]
        if '_' not in seq: 
            self.divided_chain = [(resi_list, secstruc)]
        else: 
            # cut chain into small pieces (each piece is a list of residues).
            self.divided_chain = []
            x=0
            for fr_seq, fr_secstruc in self._split_seq(seq, secstruc):
                self.divided_chain.append( (resi_list[x:x+len(fr_seq)], fr_secstruc) )
                x += len(fr_seq)
                
    def _split_seq(self, seq, secstruc):
        i_seq = 0
        i_secstruc = 0
        seq_start = 0
        secstruc_start = 0
        while i_seq < len(seq):
            if seq[i_seq] == '_':
                # generate fragment
                fr_seq = seq[seq_start:i_seq]
                fr_secstruc = secstruc[secstruc_start:i_secstruc]
                assert len(fr_seq) == len(fr_secstruc)
                yield fr_seq, fr_secstruc
                # start new fragment
                i_seq += 1
                seq_start = i_seq
                secstruc_start = i_secstruc
            i_seq += 1
            i_secstruc += 1
        # generate last fragment
        fr_seq = seq[seq_start:i_seq]
        fr_secstruc = secstruc[secstruc_start:i_secstruc]
        assert len(fr_seq) == len(fr_secstruc)
        yield fr_seq, fr_secstruc

    def divide_chain_fragment(self, chain_fragment, length):
        """
        Returns list of small lists extracted from given residues list.
        Small lists are overlapping.
        Eg. [<1>,<2>,<3>,<4>,<5>,<6>] ---> [[<1>,<2>,<3>],[<2>,<3>,<4>],[<3>,<4>,<5>][<4>,<5>,<6>]]
        
        Arguments:
        - length of small list (in the example above length = 3)
        """
        for x in range(len(chain_fragment[0])-(length-1)):
            resi_list = chain_fragment[0][x:x+length]
            secstruc = chain_fragment[1][x:x+length]
            if len(resi_list)>1: 
                yield resi_list, secstruc

    def get_fragments(self, length):
        """Generates ([resi],"secstruc") tuples."""
        for chain_fragment in self.divided_chain:
            for small_list in self.divide_chain_fragment(chain_fragment, length):
                yield small_list
            
    
class MakeLirFile:
    """
    Makes file with Lir values for all structures from given directory.

    Arguments:
    - path to file with pdb structures
    """
    def __init__(self, path=PATH_TO_LIR_STRUCTURES, list_path=PATH_TO_LIR_LIST, separator='\t'):#, output_file='LIR_DB'):
        """
        """
        self.input_data = self.get_input_list(list_path, separator)
        self.path = path
        if self.path[-1] !=os.sep: self.path+=os.sep
        self.all_records = []
        self.info = {}

    def get_input_list(self, input_file, separator='\t'):
        """Creates dict that contains names of all structures and their chains that should be included in the LIR db."""
        input_dict = {}
        f = open(input_file)
        for l in f:
            l=l.strip().split(separator)
            input_dict[l[0]] = l[1:]
        return input_dict

    def get_residue_list(self, pdb_file, chain_name):
        """
        Returns a ResidueList object.
        Arguments:
        - pdb file name
        - chain name
        """
        st = ModernaStructure('file', self.path+pdb_file, chain_name)
        return ResidueList(st)
        

    def get_lir_record(self, resi_list, pdb_filename=None,  chain_name=None, secstruc=None):
        """
        Generates record with all Lir values.
        Such record is incopleate.
        Should be compleated with structure bame and chain name.
        """
        seq = Sequence( [x.alphabet_entry for x in resi_list] )
        if not secstruc: 
            secstruc = '.'*len(seq)
        #print seq, secstruc
        l=Lir(resi_list[0], resi_list[-1])
        r = LirRecord(fr_length = len(resi_list)-2, \
                                structure = pdb_filename, \
                                chain = chain_name, \
                                preceding_resi = resi_list[0].identifier,\
                                following_resi = resi_list[-1].identifier, \
                                sequence = Sequence( seq[1:-1]) , \
                                sequence_anchor =  seq, \
                                secstruc = secstruc, \
                                x = l.x, \
                                y = l.y, \
                                dist_anchor = l.dist_anchor, \
                                beta =l.beta,  \
                                gamma = l.gamma, \
                                omega5 = l.omega5, \
                                omega3 = l.omega3, \
                                P_dist = l.P_dist, \
                                O5p_dist = l.O5p_dist, \
                                C5p_dist = l.C5p_dist,  \
                                C4p_dist = l.C4p_dist, \
                                C3p_dist = l.C3p_dist, \
                                O3p_dist = l.O3p_dist, \
                                O2p_dist = l.O2p_dist, \
                                C1p_dist = l.C1p_dist, \
                                N_dist = l.N_dist )
        return r
 
    
    def get_records_from_one_chain(self, pdb_file,  chain_name):
        """Prepares all possible LirRecords from one chain"""
        self.info[pdb_file][chain_name] = {'record_number': 0,  'record_errors': 0}
        one_chain_records_list = []
        chain_resi_list = self.get_residue_list(pdb_file,chain_name)
        for fragment_size in range(2, MAX_FRAGMENT_LENGTH):
            for small_resi_list, secstruc in chain_resi_list.get_fragments(fragment_size):
                    try:
                        new_record = self.get_lir_record(small_resi_list,  pdb_file,  chain_name, secstruc)
                        if 'x' not in new_record.sequence_anchor.seq_with_modifications:
                            one_chain_records_list.append(new_record)
                            self.info[pdb_file][chain_name] ['record_number'] +=1
                        else: self.info[pdb_file][chain_name] ['record_errors'] +=1
                    except: self.info[pdb_file][chain_name] ['record_errors'] +=1
        return one_chain_records_list


    def get_records_from_one_structure(self,  pdb_file):
        """
        Prepares all possible LirRecords for one structure.
        Returns dict: key - chain_name, value - list of LirRecords
        """
        self.info[pdb_file] = {}
        one_structure_records_list = []
        chains = self.input_data[pdb_file]#self.get_rna_chain_names(pdb_file)
        for ch in chains:
            one_structure_records_list+= self.get_records_from_one_chain(pdb_file, ch)
        return one_structure_records_list
    
    
    def generate_lir_db(self):
        """
        Generates dict with all posible records for structures from given directory
        { pdb_code : {chian_name : [residues_list] } } 
        """
        pdb_files = self.input_data.keys()#self.get_files_list()
        for pdbfile in pdb_files:
            log.write_message(pdbfile)
            self.all_records+=(self.get_records_from_one_structure(pdbfile))

    def get_lir_statistics(self):
        result = {}
        result['Number of records:'] = len(self.all_records)
        result['Number of structures:'] = len(self.info.keys())
        
        errors = 0
        chains = 0
        for pdb in self.info:
            chains +=len(self.info[pdb].keys())
            for ch in self.info[pdb]:
                errors += self.info[pdb][ch]['record_errors']
        result['Number of chains:'] = chains
        result['Number of errors:'] = errors
        
        fr_length = {}
        for x in range(0, MAX_FRAGMENT_LENGTH): fr_length[str(x)]=0
        for r in self.all_records:
           fr_length[str(r.fr_length)] +=1
        for x in range(0, MAX_FRAGMENT_LENGTH):
            result['Number of %2i residues long fragments:'%x]=fr_length[str(x)]
        return result
    

    def write_lir_db_to_file(self, db_filename='LIRdb_new_version',  log_filename='LIRdb_new_version.log'):
        """ """
        f=open(db_filename, 'w')
        #f.write ('fragment length\tstructure\tchain\tpreceding residue\tfollowing residue\tfragment sequence\tsequence with anchor\tx\ty\tdistance\tbeta\tgamma\tomega5\tomega3\n')
        f.write ('fragment length\tstructure\tchain\tpreceding residue\tfollowing residue\tfragment sequence\tsequence with anchor\tP dist\tO5p dist\tC5p dist\tC4p dist\tC3p dist\tO3p dist\tO2p dist\tC1p dist\tN* dist\n')
        for record in self.all_records:
            f.write(record.get_txt('\t')+'\n')
        f.close()
        
        g=open(log_filename, 'w')
        g.write('GENERAL FRAGMENT STATISTICS\n\n')
        statistics_dict = self.get_lir_statistics()
        statistics_dict_keys = statistics_dict.keys()
        statistics_dict_keys.sort()
        for x in statistics_dict_keys:
            g.write(x+str(statistics_dict[x])+'\n')
        g.write('\n'+50*'_'+'\n'+'DETAILED FRAGMENT STATISTICS\n\n')
        for pdb in self.info.keys():
            for ch in self.info[pdb].keys():
                g.write('Structure: %s\tChain: %s\tRecords: %s\tErrors: %s\n' %(pdb,  ch, str(self.info[pdb][ch]['record_number'] ),  str(self.info[pdb][ch]['record_errors'] ) ))
        g.close()


if __name__ == '__main__':
    log.print_all = True
    l=MakeLirFile()
    l.generate_lir_db()
    l.write_lir_db_to_file()
    
