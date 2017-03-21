#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Main lib file"""

AMINOACID_CODES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
            "TRP", "TYR", "VAL"]
RES = ['DA', 'DG', 'DT', 'DC']
RES += ['A', 'G', 'U', 'C']

RESS = ['A', 'C', 'G', 'U', 'ADE', 'CYT', 'GUA', 'URY', 'URI', 'U34', 'U31', 'C31', '4SU', 'H2U', 'QUO', 'G7M', '5MU', '5MC', 'PSU', '2MG', '1MG', '1MA', 'M2G', '5BU', 'FHU', 'FMU', 'IU', 'OMG', 'OMC', 'OMU', 'A2M', 'A23', 'CCC', 'I'] + ['RC', 'RU', 'RA', 'RG', 'RT']
#DNA = ['DA', 'DG', 'DT', 'DC']
#RNA = ['A', 'G', 'U', 'C']
IONS = ['NA', 'MG', 'MN']
HYDROGEN_NAMES = ["H", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", "H3", "H5", "H6", "H5T", "H41", "1H5'", 
                          "2H5'", "HO2'", "1H4", "2H4", "1H2", "2H2", "H1", "H8", "H2", "1H6", "2H6",
                          "HO5'", "H21", "H22", "H61", "H62", "H42", "HO3'", "1H2'", "2HO'", "HO'2", "H2'1" , "HO'2", "HO'2",
                          "H2", "H2'1", "H1", "H2", "1H5*","2H5*", "H4*", "H3*", "H1*", "1H2*", "2HO*", "1H2", "2H2", "1H4", "2H4", "1H6", "2H6", "H1", "H2", "H3", "H5", "H6", "H8", "H5'1", "H5'2", "H3T"]

import os
import sys
from collections import OrderedDict
import re
import string
import time
import urllib2
import gzip
import tempfile

from utils.extra_functions.select_fragment import select_pdb_fragment_pymol_style, select_pdb_fragment

# Don't fix OP3, ignore it
ignore_op3 = False

def get_version(currfn='', verbose=False): #dupa
    """Get version of the tool based on state of the git repository.
    Return version. 
    If currfn is empty, then the path is '.'. Hmm.. I think it will work. We will see.
    The version is not printed!
    https://github.com/m4rx9/curr_version/"""
    from commands import getoutput

    if currfn == '':
        path = '.'
    else:
        path = os.path.dirname(currfn)
    if verbose: print 'get_version::path', path
    if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
        path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))
    if not path: path = '.'
    if verbose: print 'get_version::path2', path
    curr_path = os.getcwd()
    os.chdir(os.path.abspath(path))
    version = getoutput('git describe --long --tags --dirty --always')
    os.chdir(curr_path)
    if version.find('not found')>-1:
        return ' unknown' # > install git to get versioning based on git'
    else:
        return version

class StrucFile:
    """StrucFile"""
    def __init__(self, fn):
        self.fn = fn
    
        self.report = []
        self.report.append('The RNAStrucFile report: %s ' % fn) 

        self.mol2_format = False

        self.lines = []
        lines = open(fn).read().strip().split('\n')
        has_many_models = False
        for l in lines:
            # multi-models pdb files
            if l.startswith('MODEL'):
                has_many_models = True
            if l.startswith('ENDMDL'):
                break

            if l.startswith('ATOM') or l.startswith('HETATM') or l.startswith('TER') or l.startswith('END'):

                self.lines.append(l.strip())
            if l.startswith("@<TRIPOS>"):
                self.mol2_format = True
                self.report.append('This is mol2 format')
            
        self.res = self.get_resn_uniq()

    def is_it_pdb(self):
        """Return True if the files is in PDB format."""
        if len(self.lines):
            return True
        else:
            return False

    def is_mol2(self):
        """Return True if is_mol2 based on the presence of ```@<TRIPOS>```."""
        return self.mol2_format

    def decap_gtp(self):
        lines = []
        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                if l[12:16].strip() in ['PG', 'O1G', 'O2G', 'O3G', 'O3B', 'PB','O1B','O2B', 'O3A']:
                    continue
                if l[12:16].strip() == 'PA':
                    l = l.replace('PA', 'P ')
                if l[12:16].strip() == 'O1A':
                    l = l.replace('O1A', 'O1P')
                if l[12:16].strip() == 'O2A':
                    l = l.replace('O2A', 'O2P')
                if l[17:20].strip() == 'GTP':
                    l = l[:17] + '  G' + l[20:]
                    l = l.replace('HETATM', 'ATOM  ')
                lines.append(l)
        self.lines = lines
        
    def is_amber_like(self):
        """Use self.lines and check if there is XX line
        """
        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                rn = l[17:20]
                if rn in ['RU5', 'RC5', 'RA5', 'RT5', 'RG5']:
                    self.report.append('This is amber-like format')
                    return True
        return False

    def mol2toPDB(self, outfn=""):
        try:
            import pybel
        except ImportError:
            print 'pybel is needed for mol2 to pdb convertion'
            #sys.exit(1)
            sys.exit(0)

        if not outfn:
            outfn = self.fn.replace('.mol2', '.pdb')
            
        for mol in pybel.readfile("mol2", self.fn):
            mol.write("pdb", outfn, overwrite=True)

        print 'outfn: ', outfn
        self.report.append('  Converted from mol2 to PDB')
        return outfn

    def get_no_lines(self):
        return len(self.lines)

    def get_text(self, add_end=True):
        txt = ''
        for l in self.lines:
            if l.startswith('END'):
                continue # skip end
            txt += l.strip() + '\n'
        if add_end:
            if not l.startswith('END'):
                txt += 'END'
        return txt.strip()

    def get_chain(self, chain_id='A'):
        txt = ''
        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                if l[21] == chain_id:
                    txt += l.strip() + '\n'
        txt += 'TER'
        return txt

    def get_resn_uniq(self):
        res = set()
        for l in self.lines:
            r = l[17:20].strip().upper()
            res.add(r)
        return res

    def check_res_if_std_na(self):
        wrong = []
        
        for r in self.res:
            if r not in RES:
                wrong.append(r)
        return wrong

    def get_seq(self):
        """Get seq (v2) gets segments of chains with correct numbering

        Run::

            python rna_pdb_seq.py input/1ykq_clx.pdb
            > 1ykq_clx A:101-111
            GGAGCUCGCCC
            > 1ykq_clx B:201-238
            GGGCGAGGCCGUGCCAGCUCUUCGGAGCAAUACUCGGC

            > 6_solution_0 A:1-19 26-113 117-172
            GGCGGCAGGUGCUCCCGACGUCGGGAGUUAAAAGGGAAG

        Chains is ``{'A': {'header': 'A:1-19 26-113 117-172', 'resi': [1, 2, 3, ..., \
        19, 26, 27, ..., 172], 'seq': ['G', 'G', 'C', 'G', ... C', 'G', 'U', 'C']}}``

        .. warning :: take only ATOM and HETATM lines.
        """
        seq = self.lines[0][19]
        chains = OrderedDict()
        resi_prev = None
        chain_prev = None

        chains = {}
        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                resi = int(l[22:26])
                if resi_prev != resi:
                    resname = l[17:20].strip()
                    chain_curr = l[21]

                    if len(resname) == 'GTP': # DG -> g GTP
                        resname = 'g'
                    if len(resname) > 1: # DG -> g GTP
                        resname = resname[-1].lower()

                    try:
                        chains[chain_curr]['resi'].append(resi)
                        chains[chain_curr]['seq'].append(resname)
                    except KeyError:
                        chains[chain_curr] = {}
                        chains[chain_curr]['resi'] = []
                        chains[chain_curr]['resi'].append(resi)
                        chains[chain_curr]['seq'] = []
                        chains[chain_curr]['seq'].append(resname)

                    resi_prev = resi
                    chain_prev = chain_curr

        for c in chains.keys():
            header = c + ':' + str(chains[c]['resi'][0]) + '-' # add A:1-
            for i in range(1, len(chains[c]['resi'])): # start from second element
                if chains[c]['resi'][i] - chains[c]['resi'][i - 1] > 1:
                    header += '' + str(chains[c]['resi'][i - 1]) + ' '
                    header += '' + str(chains[c]['resi'][i]) + '-'                
            header += '' + str(chains[c]['resi'][-1])                
            chains[c]['header'] = header # add -163 (last element)

        txt = ''
        for c in chains.keys():
            txt += '> ' + os.path.basename(self.fn.replace('.pdb', '')) + ' ' + chains[c]['header'] + '\n'
            txt += ''.join(chains[c]['seq']) + '\n'
        return txt.strip()


    def __get_seq(self):
        """get_seq DEPRECATED

        You get `chains` such as:
        OrderedDict([('A', {'header': 'A:1-47', 'seq': 'CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAU'}), ('B', {'header': 'B:48-58', 'seq': 'AUCAGGUGCAA'})])

        .. warning:: take only ATOM and HETATM lines.
        """
        seq = ''
        curri = int(self.lines[0][22:26])
        seq = self.lines[0][19]
        chains = OrderedDict()
        curri = -100000000000000000000000000000000 #ugly
        chain_prev = None

        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                resi = int(l[22:26])
                if curri != resi:
                    print l
                    resname = l[17:20].strip()
                    if len(resname) == 'GTP': # DG -> g GTP
                        resname = 'g'
                    if len(resname) > 1: # DG -> g GTP
                        resname = resname[-1].lower()

                    seq += resname
                    chain_curr = l[21]

                    # if distances between curr res and previevs is bigger than 1, then show it as a fragment
                    if resi - curri > 1 and resi - curri < 100000000000000000000000000000000: # ugly hack
                        print resi - curri
                        chains[chain_prev]['header'] += '-' + str(resi_prev)                        
                    if chain_prev != chain_curr and chain_prev:
                        chains[chain_prev]['header'] += '-' + str(resi_prev)
                    if chains.has_key(chaoin_curr):                
                        chains[chain_curr]['seq'] += resname
                    else:
                        chains[chain_curr] = dict()
                        chains[chain_curr]['header'] = chain_curr + ':' + str(resi)#resi_prev)
                        chains[chain_curr]['seq'] = resname
                    resi_prev = resi
                    chain_prev = chain_curr
                curri = resi
        chains[chain_prev]['header'] += '-' + str(resi_prev)
        seq = ''
        for c in chains.keys():
            seq += '> ' + os.path.basename(self.fn) + ' ' + chains[c]['header'] + '\n'
            seq += chains[c]['seq'] + '\n'
        return seq.strip()

    def get_info_chains(self):
        """return A:3-21 B:22-32
        """
        seq = ''
        curri = int(self.lines[0][22:26])
        seq = self.lines[0][19]
        chains = OrderedDict()
        curri = -100000000000000 #ugly
        chain_prev = None
        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                resi = int(l[22:26])
                if curri != resi:
                    resname = l[17:20].strip()
                    if len(resname) == 'GTP': # DG -> g GTP
                        resname = 'g'
                    if len(resname) > 1: # DG -> g GTP
                        resname = resname[-1].lower()
                    seq += resname
                    chain_curr = l[21]
                    if chain_prev != chain_curr and chain_prev:
                        chains[chain_prev]['header'] += '-' + str(resi_prev)
                    if chains.has_key(chain_curr):                
                        chains[chain_curr]['seq'] += resname
                    else:
                        chains[chain_curr] = dict()
                        chains[chain_curr]['header'] = chain_curr + ':' + str(resi)#resi_prev)
                        chains[chain_curr]['seq'] = resname
                    resi_prev = resi
                    chain_prev = chain_curr
                curri = resi
        chains[chain_prev]['header'] += '-' + str(resi_prev)
        seq = ''
        for c in chains.keys():
            seq += chains[c]['header'] + ' '
        return seq.strip()

    def detect_file_format(self):
        pass
    
    def detect_molecule_type(self):
        aa = []
        na = []
        for r in self.res:
            if r in AMINOACID_CODES:
                aa.append(r)
            if r in RESS:
                na.append(r)            

        aa = float(len(aa)) / len(self.res)
        na = float(len(na)) / len(self.res)

        if aa == 0 and na == 0:
            return 'error'
        if aa > na:
            return '>protein< vs na', aa, na
        else:
            return 'protein vs >na<', aa, na

    def get_head(self):
        return '\n'.join(self.lines[:5])

    def get_tail(self):
        return '\n'.join(self.lines[-5:])

    def get_preview(self):
        t = '\n'.join(self.lines[:5])
        t += '\n-------------------------------------------------------------------\n'
        t += '\n'.join(self.lines[-5:])
        return t

    def remove_hydrogen(self):
        lines = []
        for l in self.lines:
            if l[77:79].strip() == 'H':
                continue
            if l[12:16].strip() in HYDROGEN_NAMES:
            #if l[12:16].strip().startswith('H'):
                continue
            else:
                #print l[12:16]
                lines.append(l)
        self.lines = lines

    def remove_water(self):
        """Remove HOH and TIP3"""
        lines  = []
        for l in self.lines:
            if l[17:21].strip() in ['HOH', 'TIP3', 'WAT']:
                continue
            else:
                lines.append(l)
        self.lines = lines

    def remove_ion(self):
        """
    TER    1025        U A  47                                                      
    HETATM 1026 MG    MG A 101      42.664  34.395  50.249  1.00 70.99          MG  
    HETATM 1027 MG    MG A 201      47.865  33.919  48.090  1.00 67.09          MG 
        """
        lines = []
        for l in self.lines:
            element = l[76:78].strip().upper()
            element2 = l[17:20].strip().upper()
            if element in IONS:
                continue
            if element2 in IONS:
                continue
            else:
                lines.append(l)
        self.lines = lines

    def fixU__to__U(self):
        lines = []
        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                rn = l[17:20]
                rn = rn.replace('G  ', '  G')
                rn = rn.replace('U  ', '  U')
                rn = rn.replace('C  ', '  C')
                rn = rn.replace('A  ', '  A')
                l = l[:16] + ' ' + rn + ' ' + l[21:]
            #print l.strip()
            #print l2
            #l = l.replace(' U   ', '   U ')
            #l = l.replace(' G   ', '   G ')
            #l = l.replace(' A   ', '   A ')
            #l = l.replace(' C   ', '   C ')
            lines.append(l)
        print 'fixU__to__U OK'
        self.report.append('  Fix: U__ -> __U')
        self.lines = lines

    def resn_as_dna(self):
        lines = []
        for l in self.lines:
            if l.startswith('ATOM') or l.startswith('HETATM') :
                #print l
                nl = l.replace( 'DA5', ' DA') # RA should be the last!!!!
                nl = nl.replace('DA3', ' DA')
                nl = nl.replace(' DA', ' DA')
                nl = nl.replace(' rA', ' DA')

                nl = nl.replace('DC5', ' DC')
                nl = nl.replace('DC3', ' DC')
                nl = nl.replace(' DC', ' DC')
                nl = nl.replace(' rC', ' DC')

                nl = nl.replace('DG5', ' DG')
                nl = nl.replace('DG3', ' DG')
                nl = nl.replace(' DG', ' DG')
                nl = nl.replace(' rG', ' DG')

                nl = nl.replace('DU5', ' DU')
                nl = nl.replace('DU3', ' DU')
                nl = nl.replace(' DU', ' DU')
                nl = nl.replace(' rU', ' DU')

                nl = nl.replace('DT5', ' DT')
                nl = nl.replace('DT3', ' DT')
                nl = nl.replace(' DT', ' DT')
                nl = nl.replace(' rT', ' DT')

                nl = nl.replace('C5M', 'C7 ')

                if l[17:20].strip() == 'G':
                    nl = nl[:17] + ' DG' + nl[20:]

                if l[17:20].strip() == 'C':
                    nl = nl[:17] + ' DC' + nl[20:]

                if l[17:20].strip() == 'T':
                    nl = nl[:17] + ' DT' + nl[20:]

                if l[17:20].strip() == 'U':
                    nl = nl[:17] + ' DU' + nl[20:]

                if l[17:20].strip() == 'A':
                    nl = nl[:17] + ' DA' + nl[20:]

                lines.append(nl)
            if l.startswith("END") or l.startswith("TER"):
                lines.append(l)

        print 'resn_as_dna'
        self.report.append('  resn_as_dna')
        self.lines = lines

    def fix_O_in_UC(self):
        """.. warning: remove RU names before using this function"""
        lines = []
        for l in self.lines:
            #if l[12:16].strip() in 
            #if l[12:16].strip().startswith('H'):
            nl = l.replace('O     U',
                           'O2    U')
            nl =nl.replace('O     C',
                           'O2    C')
            lines.append(nl)
        self.lines = lines


    def fix_op_atoms(self):
        """Replace OXP' to OPX1, e.g ('O1P' -> 'OP1')"""
        lines = []
        for l in self.lines:
            nl = l.replace('*', '\'')
            nl = nl.replace('O1P', 'OP1')
            nl = nl.replace('O2P', 'OP2')
            nl = nl.replace('O3P', 'OP3')
            lines.append(nl)
        self.lines = lines

    def get_report(self):
        """
        :return report: string
        """
        return '\n'.join(self.report)

    def is_rna(self):
        wrong = []
        for r in self.res:
            if r.upper().strip() in ['RC', 'RU', 'RA', 'RG', 'RT']:
                if r not in wrong_res:
                    wrong_res.append(r)
        return wrong_res

    def check_res_if_std_dna(self):
        wrong_res = []
        for r in self.res:
            if r.upper().strip() in ['A', 'T', 'C', 'G']:
                if r not in wrong_res:
                    wrong_res.append(r)
        return wrong_res

    def check_res_if_supid_rna(self):
        wrong_res = []
        for r in self.res:
            if r.upper().strip() in ['RC', 'RU', 'RA', 'RG', 'RT']:
                if r not in wrong_res:
                    wrong_res.append(r)
        return wrong_res

    def is_rna(self):
        for r in self.res:
            if r.upper().strip() in ['RC', 'RU', 'RA', 'RG', 'RT']:
                if r not in wrong_res:
                    wrong_res.append(r)
        return wrong_res

    def renum_atoms(self):
        """Renum atoms, from 1 to X for line; ATOM/HETATM"""
        lines = []
        c = 1
        for l in self.lines:
            l = l[:6] + str(c).rjust(5) + l[11:]
            c += 1
            lines.append(l)
        self.lines = lines

    def fix_resn(self):
        """
        fix::

         # URI -> U, URA -> U
         1xjr_clx_charmm.pdb:ATOM    101  P   URA A   5      58.180  39.153  30.336  1.00 70.94
         rp13_Dokholyan_1_URI_CYT_ADE_GUA_hydrogens.pdb:ATOM  82  P   URI A   4     501.633 506.561 506.256  1.00  0.00           P"""
        lines = []
        for l in self.lines:
            nl = l.replace( 'RA5', '  A') # RA should be the last!!!!
            nl = nl.replace('RA3', '  A')
            nl = nl.replace('ADE', '  A')
            nl = nl.replace(' RA', '  A')
            nl = nl.replace(' rA', '  A')

            nl = nl.replace('RC5', '  C')
            nl = nl.replace('RC3', '  C')
            nl = nl.replace('CYT', '  C')
            nl = nl.replace(' RC', '  C')
            nl = nl.replace(' rC', '  C')

            nl = nl.replace('RG5', '  G')
            nl = nl.replace('RG3', '  G')
            nl = nl.replace('GUA', '  G')
            nl = nl.replace(' RG', '  G')
            nl = nl.replace(' rG', '  G')

            nl = nl.replace('RU5', '  U')
            nl = nl.replace('RU3', '  U')
            nl = nl.replace('URA', '  U')
            nl = nl.replace('URI', '  U')              
            nl = nl.replace('URY', '  U')              
            nl = nl.replace(' RU', '  U')
            nl = nl.replace(' rU', '  U')

            nl = nl.replace('RT5', '  T')
            nl = nl.replace('RT3', '  T')
            nl = nl.replace('THY', '  T')  
            nl = nl.replace(' RT', '  T')
            nl = nl.replace(' rT', '  T')
            
            lines.append(nl)
            
        self.lines = lines

    def check_res_if_std_prot(self):
        wrong = []
        for r in self.res:
            if r not in AMINOACID_CODES:
                wrong.append(r)
        return wrong

    def write(self, outfn,v=True):
        """Write ```self.lines``` to a file (and END file")"""
        f = open(outfn, 'w')
        for l in self.lines:
            f.write(l + '\n')
        if not l.startswith('END'):
            f.write('END')
        f.close()
        if v:
            print 'Write %s' % outfn

    def get_atom_num(self,line):
        """Extract atom number from a line of PDB file
        Arguments:
          * line = ATOM line from a PDB file
        Output:
          * atom number as an integer
        """
        return int(''.join(filter(lambda x: x.isdigit(), line[6:11]))) 

    def get_res_num(self,line):
        """Extract residue number from a line of PDB file
        Arguments:
          * line = ATOM line from a PDB file
        Output:
          * residue number as an integer
        """
        return int(''.join(filter(lambda x: x.isdigit(), line[22:27])))

    def get_res_code(self,line):
        """Get residue code from a line of a PDB file
        """
        if not line.startswith('ATOM'):
            return None
        return line[17:20]

    def get_atom_code(self,line):
        """Get atom code from a line of a PDB file
        """
        if not line.startswith('ATOM'):
            return None
        return line[13:16].replace(' ', '')

    def get_atom_coords(self,line):
        """Get atom coordinates from a line of a PDB file
        """
        if not line.startswith('ATOM'):
            return None
        return tuple(map(float, line[31:54].split()))

    def set_line_bfactor(self,line, bfactor):
        if not line.startswith('ATOM'):
            return None
        return line[:60] + (" %5.2f" % bfactor) + line[66:]

    def set_atom_occupancy(self, line, occupancy):
        """set occupancy for line"""
        return line[:54] + (" %5.2f" % occupancy) + line[60:]

    def set_atom_code(self,line, code):
        return line[:13] + code + ' ' * (3 - len(code)) + line[16:]

    def set_res_code(self,line, code):
        return line[:17] + code.rjust(3) + line[21:]

    def get_chain_id(self,line):
        return line[21:22]

    def get_atom_index(self,line):
        try:
            return int(line[6:11])
        except:
            return None

    def set_atom_index(self,line,index):
        return line[:7] + str(index).rjust(4) + line[11:]

    def get_res_index(self,line):
        return int(line[22:26])

    def set_res_index(self,line, index):
        return line[:23] + str(index).rjust(3) + line[26:]

    def set_chain_id(self,line, chain_id):
        return line[:21] + chain_id + line[22:]

    def get_rnapuzzle_ready(self, renumber_residues=True, fix_missing_atoms=False):#:, ready_for="RNAPuzzle"):
        """Get rnapuzzle (SimRNA) ready structure.

        Clean up a structure, get corrent order of atoms.

        :param renumber_residues: boolean, from 1 to ..., second chain starts from 1 etc. 
        :param fix_missing_atoms: boolean, superimpose motifs from the minilibrary and copy-paste missing atoms, this is super crude, so should be used with caution.

        Submission format @http://ahsoka.u-strasbg.fr/rnapuzzles/

        Run :func:`rna_pdb_tools.pdb_parser_lib.StrucFile.fix_resn` before this function to fix names.

        - 170305 Merged with get_simrna_ready and fixing OP3 terminal added 
        - 170308 Fix missing atoms for bases, and O2'

        .. image:: ../pngs/fix_missing_o_before_after.png
        Fig. Add missing O2' atom (before and after).

        .. image:: ../pngs/fix_missing_superposition.png
        Fig. The residue to fix is in cyan. The G base from the library in red. Atoms O4', C2', C1' are shared between the sugar (in cyan) and the G base from the library (in red). These atoms are used to superimpose the G base on the sugar, and then all atoms from the base are copied to the residues.

        .. image:: ../pngs/fix_missing_bases.png
        **Fig.** Rebuild ACGU base-less. It's not perfect but good enough for some applications.

        .. warning:: It was only tested with the whole base missing!
        
        .. warning:: requires: Biopython"""
        try:
            from Bio import PDB
            from Bio.PDB import PDBIO
            import warnings
            warnings.filterwarnings('ignore', '.*Invalid or missing.*',)
            warnings.filterwarnings('ignore', '.*with given element *',)
        except:
            sys.exit('Error: Install biopython to use this function (pip biopython)')

        import copy
        # for debugging
        #v = True
        v = False
        #renumber_residues = True
        
        #if ready_for == "RNAPuzzle":
        G_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split()
        A_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split()
        U_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split()
        C_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split()

        # hmm.. is it the same as RNApuzzle???
        #if ready_for == "SimRNA":
        #    G_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split()
        #    A_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split()
        #    U_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split()
        #    C_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split()

        tf = tempfile.NamedTemporaryFile(delete=False)
        ftmp = tf.name
        self.write(ftmp, v=False)

        parser = PDB.PDBParser()
        struct = parser.get_structure('', ftmp)
        model = struct[0]

        s2 = PDB.Structure.Structure(struct.id)
        m2 = PDB.Model.Model(model.id)

        chains2 = []

        missing = []
        fixed = []
        
        for chain in model.get_list():
            if v: print 'chain:', chain
            res = [] 
            for r in chain:
                res.append(r)

            res = copy.copy(res)

            c2 = PDB.Chain.Chain(chain.id)        

            c = 1  # new chain, goes from 1 !!! if renumber True
            for r in res:
                # hack for amber/qrna
                r.resname = r.resname.strip()
                if r.resname == 'RC3': r.resname = 'C'
                if r.resname == 'RU3': r.resname = 'U'
                if r.resname == 'RG3': r.resname = 'G'
                if r.resname == 'RA3': r.resname = 'A'

                if r.resname == 'C3': r.resname = 'C'
                if r.resname == 'U3': r.resname = 'U'
                if r.resname == 'G3': r.resname = 'G'
                if r.resname == 'A3': r.resname = 'A'

                if r.resname == 'RC5': r.resname = 'C'
                if r.resname == 'RU5': r.resname = 'U'
                if r.resname == 'RG5': r.resname = 'G'
                if r.resname == 'RA5': r.resname = 'A'

                if r.resname == 'C5': r.resname = 'C'
                if r.resname == 'U5': r.resname = 'U'
                if r.resname == 'G5': r.resname = 'G'
                if r.resname == 'A5': r.resname = 'A'

                if r.resname.strip() == 'RC': r.resname = 'C'
                if r.resname.strip() == 'RU': r.resname = 'U'
                if r.resname.strip() == 'RG': r.resname = 'G'
                if r.resname.strip() == 'RA': r.resname = 'A'

                r2 = PDB.Residue.Residue(r.id, r.resname.strip(), r.segid)
                if renumber_residues:
                    r2.id = (r2.id[0], c, r2.id[2]) ## renumber residues
                #
                # experimental: fixing missing OP3
                #
                if c == 1:
                    # if p_missing
                    p_missing = True
                    #if p_missing:
                    #    try:
                    #        x = r["O5'"]
                    #        x.id =       ' P'
                    #        x.name =     ' P'
                    #        x.fullname = ' P'
                    #        print "REMARK 000 FIX O5' -> P fix in chain ", chain.id
                    #    except:
                    #        pass
                    for a in r:
                        if a.id == 'P':
                            p_missing = False
                    if v: print 'p_missing', p_missing

                    if p_missing and fix_missing_atoms:
                            currfn = __file__
                            if currfn == '':
                                path = '.'
                            else:
                                path = os.path.dirname(currfn)
                            if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
                                path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))

                            po3_struc = PDB.PDBParser().get_structure('', path + '/data/PO3_inner.pdb') 
                            po3 = [po3_atom for po3_atom in po3_struc[0].get_residues()][0]

                            r_atoms = [r["O4'"], r["C4'"], r["C3'"]]
                            po3_atoms = [po3["O4'"], po3["C4'"], po3["C3'"]]

                            sup = PDB.Superimposer()
                            sup.set_atoms(r_atoms, po3_atoms)
                            rms = round(sup.rms, 3)

                            sup.apply( po3_struc.get_atoms() ) # to all atoms of po3

                            r.add( po3['P'])
                            r.add( po3['OP1'])
                            r.add( po3['OP2'])
                            try:
                                r.add( po3["O5'"]) 
                            except:
                                del r["O5'"] 
                                r.add( po3["O5'"]) 

                    p_missing = False # off this function

                    # save it
                    #io = PDB.PDBIO()
                    #io.set_structure( po3_struc )
                    #io.save("po3.pdb")

                    # if o2p_missing
                    o2p_missing = True
                    #if p_missing:
                    for a in r:
                        if a.id == "O2'":
                            o2p_missing = False
                    if v: print 'o2p_missing', o2p_missing

                    if o2p_missing and fix_missing_atoms:
                            currfn = __file__
                            if currfn == '':
                                path = '.'
                            else:
                                path = os.path.dirname(currfn)
                            if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
                                path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))

                            o2p_struc = PDB.PDBParser().get_structure('', path + '/data/o2prim.pdb') 
                            o2p = [o2p_atom for o2p_atom in o2p_struc[0].get_residues()][0]

                            r_atoms = [r["C3'"], r["C2'"], r["C1'"]]
                            o2p_atoms = [o2p["C3'"], o2p["C2'"], o2p["C1'"]]

                            sup = PDB.Superimposer()
                            sup.set_atoms(r_atoms, o2p_atoms)
                            rms = round(sup.rms, 3)

                            sup.apply( o2p_struc.get_atoms() ) # to all atoms of o2p

                            r.add( o2p["O2'"])

                    o2p_missing = False # off this function

                    # fix C
                if str(r.get_resname()).strip() == "C" and fix_missing_atoms:
                    for a in r:
                        if a.id == "N1":
                            break
                    else: # fix
                            currfn = __file__
                            if currfn == '':
                                path = '.'
                            else:
                                path = os.path.dirname(currfn)
                            if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
                                path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))

                            C_struc = PDB.PDBParser().get_structure('', path + '/data/C.pdb') 
                            C = [C_atom for C_atom in C_struc[0].get_residues()][0]

                            r_atoms = [r["O4'"], r["C2'"], r["C1'"]]
                            C_atoms = [C["O4'"], C["C2'"], C["C1'"]]

                            sup = PDB.Superimposer()
                            sup.set_atoms(r_atoms, C_atoms)
                            rms = round(sup.rms, 3)

                            sup.apply( C_struc.get_atoms() ) # to all atoms of C

                            r.add( C["N1"])
                            r.add( C["C2"])
                            r.add( C["O2"])
                            r.add( C["N3"])
                            r.add( C["C4"])                            
                            r.add( C["N4"])
                            r.add( C["C5"])
                            r.add( C["C6"])

                            fixed.append(['add the whole base C', chain.id, r, c])
                            
                # fix U
                if str(r.get_resname()).strip() == "U" and fix_missing_atoms:
                    for a in r:
                        if a.id == "N1":
                            break
                    else: # fix
                            currfn = __file__
                            if currfn == '':
                                path = '.'
                            else:
                                path = os.path.dirname(currfn)
                            if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
                                path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))

                            U_struc = PDB.PDBParser().get_structure('', path + '/data/U.pdb') 
                            U = [U_atom for U_atom in U_struc[0].get_residues()][0]

                            r_atoms = [r["O4'"], r["C2'"], r["C1'"]]
                            U_atoms = [U["O4'"], U["C2'"], U["C1'"]]

                            sup = PDB.Superimposer()
                            sup.set_atoms(r_atoms, U_atoms)
                            rms = round(sup.rms, 3)

                            sup.apply( U_struc.get_atoms() ) # to all atoms of U

                            r.add( U["N1"])
                            r.add( U["C2"])
                            r.add( U["O2"])
                            r.add( U["N3"])
                            r.add( U["C4"])                            
                            r.add( U["O4"])
                            r.add( U["C5"])
                            r.add( U["C6"])

                            fixed.append(['add the whole base U', chain.id, r, c])
                # fix G
                if str(r.get_resname()).strip() == "G" and fix_missing_atoms:
                    for a in r:
                        if a.id == "N1":
                            break
                    else: # fix
                            currfn = __file__
                            if currfn == '':
                                path = '.'
                            else:
                                path = os.path.dirname(currfn)
                            if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
                                path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))

                            G_struc = PDB.PDBParser().get_structure('', path + '/data/G.pdb') 
                            G = [G_atom for G_atom in G_struc[0].get_residues()][0]

                            r_atoms = [r["O4'"], r["C2'"], r["C1'"]]
                            G_atoms = [G["O4'"], G["C2'"], G["C1'"]]

                            sup = PDB.Superimposer()
                            sup.set_atoms(r_atoms, G_atoms)
                            rms = round(sup.rms, 3)

                            sup.apply( G_struc.get_atoms() ) # to all atoms of G

                            r.add( G["N9"])
                            r.add( G["C8"])
                            r.add( G["N7"])
                            r.add( G["C5"])
                            r.add( G["C6"])                            
                            r.add( G["O6"])
                            r.add( G["N1"])
                            r.add( G["C2"])
                            r.add( G["N2"])
                            r.add( G["N3"])
                            r.add( G["C4"])

                            fixed.append(['add the whole base G', chain.id, r, c])
                # fix A
                if str(r.get_resname()).strip() == "A" and fix_missing_atoms:
                    for a in r:
                        if a.id == "N1":
                            break
                    else: # fix
                            currfn = __file__
                            if currfn == '':
                                path = '.'
                            else:
                                path = os.path.dirname(currfn)
                            if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
                                path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))

                            A_struc = PDB.PDBParser().get_structure('', path + '/data/A.pdb') 
                            A = [A_atom for A_atom in A_struc[0].get_residues()][0]

                            r_atoms = [r["O4'"], r["C2'"], r["C1'"]]
                            A_atoms = [A["O4'"], A["C2'"], A["C1'"]]

                            sup = PDB.Superimposer()
                            sup.set_atoms(r_atoms, A_atoms)
                            rms = round(sup.rms, 3)

                            sup.apply( A_struc.get_atoms() ) # to all atoms of A

                            r.add( A["N9"])
                            r.add( A["C8"])
                            r.add( A["N7"])
                            r.add( A["C5"])
                            r.add( A["C6"])                            
                            r.add( A["N6"])
                            r.add( A["N1"])
                            r.add( A["C2"])
                            r.add( A["N3"])
                            r.add( A["C4"])

                            fixed.append(['add the whole base A', chain.id, r, c])

                if str(r.get_resname()).strip() == "G":
                    for an in G_ATOMS:
                        if c == 1 and ignore_op3:
                            if an in ['P', 'OP1', 'OP2']:
                                continue
                        try:
                            if c == 1 and an == "O5'" and p_missing:
                                r2.add(x)
                            else:
                                r2.add(r[an])
                        except KeyError:
                            #print 'Missing:', an, r, ' new resi', c
                            missing.append([an, chain.id, r, c])
                    c2.add(r2)

                elif str(r.get_resname()).strip() == "A":
                    for an in A_ATOMS:
                        if c == 1 and ignore_op3:
                            if an in ['P', 'OP1', 'OP2']:
                                continue
                        try:
                            if c == 1 and an == "O5'" and p_missing:
                                r2.add(x)
                            else:
                                r2.add(r[an])
                        except KeyError:
                            #print 'Missing:', an, r, ' new resi', c
                            missing.append([an, chain.id, r, c])
                    c2.add(r2)

                elif str(r.get_resname()).strip() == "C":
                    for an in C_ATOMS:
                        if c == 1 and ignore_op3:
                            if an in ['P', 'OP1', 'OP2']:
                                continue
                        try:
                            if c == 1 and an == "O5'" and p_missing:
                                r2.add(x)
                            else:
                                r2.add(r[an])
                        except:
                            #print 'Missing:', an, r, ' new resi', c
                            missing.append([an, chain.id, r, c])
                    c2.add(r2)

                elif str(r.get_resname()).strip() == "U":
                    for an in U_ATOMS:
                        if c == 1 and ignore_op3:
                            if an in ['P', 'OP1', 'OP2']:
                                continue
                        try:
                            if c == 1 and an == "O5'" and p_missing:
                                r2.add(x)
                            else:
                                r2.add(r[an])
                        except KeyError:
                            #print 'Missing:', an, r,' new resi', c
                            missing.append([an, chain.id, r, c])
                    c2.add(r2)

                c += 1
            chains2.append(c2)

        io = PDBIO()
        s2.add(m2)
        for chain2 in chains2:
            m2.add(chain2) 
        #print c2
        #print m2
        io.set_structure(s2)

        tf = tempfile.NamedTemporaryFile(delete=False)
        fout = tf.name
        io.save(fout)
        
        if fixed:
            print 'REMARK 000 Fixed atoms/residues:'
            for i in fixed:
                print 'REMARK 000 - ', i[0], 'in chain:', i[1], i[2], 'residue #', i[3]
                
        if missing:
            print 'REMARK 000 Missing atoms:'
            for i in missing:
                print 'REMARK 000  +', i[0], i[1], i[2], 'residue #', i[3]
            #raise Exception('Missing atoms in %s' % self.fn)
        #
        # fix ter 'TER' -> TER    1528        G A  71
        #
        s = StrucFile(fout)
        self.lines = s.lines
        c = 0
        #ATOM   1527  C4    G A  71       0.000   0.000   0.000  1.00  0.00           C
        nlines = []
        no_ters = 0
        for l in self.lines:
            if l.startswith('TER'):
                atom_l = self.lines[c-1]
                #print 'TER    1528        G A  71 <<<'
                new_l = 'TER'.ljust(80)
                new_l = self.set_atom_index(new_l, str(self.get_atom_index(atom_l)+1 + no_ters))
                new_l = self.set_res_code(new_l, self.get_res_code(atom_l))
                new_l = self.set_chain_id(new_l, self.get_chain_id(atom_l))
                new_l = self.set_res_index(new_l, self.get_res_index(atom_l))                
                #print new_l
                nlines.append(new_l)
                no_ters += 1
            else:
                if self.get_atom_index(l):
                    l = self.set_atom_index(l, self.get_atom_index(l) + no_ters) # 1 ter +1 2 ters +2 etc
                nlines.append(l)
            c += 1
        self.lines = nlines

    def set_occupancy_atoms(self, occupancy):
        """
        :param occupancy:
        """
        nlines = []
        for l in self.lines:
           if l.startswith('ATOM'):
               l = self.set_atom_occupancy(l, 0.00)
               nlines.append(l)
           else:
               nlines.append(l)
        self.lines = nlines
               
    def edit_occupancy_of_pdb(txt, pdb, pdb_out,v=False):
        """Make all atoms 1 (flexi) and then set occupancy 0 for seletected atoms.
        Return False if error. True if OK
        """
        struc = PDB.PDBParser().get_structure('struc', pdb)

        txt = txt.replace(' ','')
        if v:print txt
        l = re.split('[,:;]', txt)
        if v:print l 

        for s in struc:
            for c in s:
                for r in c:
                    for a in r:
                        a.set_occupancy(1)  # make it flaxi

        for i in l: # ['A', '1-10', '15', '25-30', 'B', '1-10']

            if i in string.ascii_letters:
                if v:print 'chain', i
                chain_curr = i
                continue

            if i.find('-') > -1:
                start, ends = i.split('-')
                if start > ends:
                    print >>sys.stderr, 'Error: range start > end ' + i
                    return False
                index = range(int(start), int(ends)+1)
            else:
                index=[int(i)]

            for i in index:
                # change b_factor
                try:
                    atoms = struc[0][chain_curr][i]
                except KeyError:
                    if i == chain_curr:
                        print >>sys.stderr, 'Error: Chain ' + chain_curr + ' not found in the PDB structure'
                    else:
                        print >>sys.stderr, 'Error: Residue ' + chain_curr + ':' + str(i) + ' found in the PDB structure'
                        return False
                for a in atoms:
                    a.set_occupancy(0)

        io = PDBIO()
        io.set_structure(struc)
        io.save(pdb_out)
        print 'Saved ', pdb_out
        return True

    def view(self):
        os.system('pymol ' + self.fn)


def add_header(version=None):
    now = time.strftime("%c")
    print 'HEADER Generated with rna-pdb-tools'
    print 'HEADER ver %s \nHEADER https://github.com/mmagnus/rna-pdb-tools \nHEADER %s' % (version, now)

def edit_pdb(args):
    """Edit your structure. 

    The function can take ``A:3-21>A:1-19`` or even syntax like this
    ``A:3-21>A:1-19,B:22-32>B:20-30`` and will do an editing.

    The output is printed, line by line. Only ATOM lines are edited!

    Examples::

      $ rna_pdb_tools.py --edit 'A:3-21>A:1-19' 1f27_clean.pdb > 1f27_clean_A1-19.pdb

    or even::

      $ rna_pdb_tools.py --edit 'A:3-21>A:1-19,B:22-32>B:20-30' 1f27_clean.pdb > 1f27_clean_renumb.pdb

    """
    ## open a new file
    s = StrucFile(args.file)
    if not args.no_hr:
        add_header()
        print 'HEADER --edit ' + args.edit
        
    ## --edit 'A:3-21>A:1-19,B:22-32>B:20-30'
    if args.edit.find(',')>-1:
        # more than one edits
        edits = args.edit.split(',') # ['A:3-21>A:1-19', 'B:22-32>B:20-30']
        selects = []
        for e in edits:
            selection_from, selection_to = select_pdb_fragment(e.split('>')[0]), select_pdb_fragment(e.split('>')[1])
            if len(selection_to) != len(selection_from):
                raise Exception('len(selection_to) != len(selection_from)')
            selects.append([selection_from, selection_to])
        print edits
    else:
        # one edit
        e = args.edit
        selection_from, selection_to = select_pdb_fragment(e.split('>')[0]), select_pdb_fragment(e.split('>')[1])
        if len(selection_to) != len(selection_from):
            raise Exception('len(selection_to) != len(selection_from)')
        selects = [[selection_from, selection_to]]

    ## go ever all edits: ['A:3-21>A:1-19','B:22-32>B:20-30']
    for l in s.lines:
            if l.startswith('ATOM'):
                # get chain and resi
                chain = l[21:22].strip()
                resi = int(l[22:26].strip())

                if_selected_dont_print = False
                # for selections
                for select in selects:
                    selection_from, selection_to = select
                    if selection_from.has_key(chain):
                        if resi in selection_from[chain]:
                            # [1,2,3] mapping from [4,5,10], you want to know how to map 1
                            # 1 is [0] element of first list, so you have to index first list
                            # to get 0, with this 0 you can get 4 out of second list [4,5,10][0] -> 4
                            nl = list(l)
                            chain_new = selection_to.keys()[0] # chain form second list
                            nl[21] =  chain_new # new chain
                            index = selection_from[chain].index(int(resi)) # get index of 1
                            resi_new = str(selection_to[chain_new][index]).rjust(4) # 'A' [1,2,3] -> '  1'
                            nl[22:26] = resi_new
                            nl = ''.join(nl)
                            if_selected_dont_print = True
                            print nl
                if not if_selected_dont_print:
                    print l
            else: # if not atom
                print l
    
def collapsed_view(args):
    """Collapsed view of pdb file. Only lines with C5' atoms are shown and TER, MODEL, END.
    
    example::

        [mm] rna_pdb_tools git:(master) $ python rna-pdb-tools.py --cv input/1f27.pdb
        ATOM      1  C5'   A A   3      25.674  19.091   3.459  1.00 16.99           C
        ATOM     23  C5'   C A   4      19.700  19.206   5.034  1.00 12.65           C
        ATOM     43  C5'   C A   5      14.537  16.130   6.444  1.00  8.74           C
        ATOM     63  C5'   G A   6      11.726  11.579   9.544  1.00  9.81           C
        ATOM     86  C5'   U A   7      12.007   7.281  13.726  1.00 11.35           C
        ATOM    106  C5'   C A   8      12.087   6.601  18.999  1.00 12.74           C
        TER""" 
    r = StrucFile(args.file)
    for l in r.lines:
        at = r.get_atom_code(l)
        if  at == "C5'":
            print l
        if l.startswith('TER') or l.startswith('MODEL') or l.startswith('END'):
            print l

def fetch(pdb_id):
    """fetch pdb file from RCSB.org
    https://files.rcsb.org/download/1Y26.pdb"""
    try:
        response = urllib2.urlopen('https://files.rcsb.org/download/' + pdb_id + '.pdb')
    except urllib2.HTTPError:
        raise Exception('The PDB does not exists: ' + pdb_id)
    txt = response.read()

    print 'downloading...' + pdb_id,
    
    with open(pdb_id + '.pdb', 'w') as f:
        f.write(txt)
    print 'ok'

def fetch_ba(pdb_id):
    """fetch biological assembly pdb file from RCSB.org"""
    try:
        response = urllib2.urlopen('https://files.rcsb.org/download/' + pdb_id + '.pdb1.gz')
    except urllib2.HTTPError:
        raise Exception('The PDB does not exists: ' + pdb_id)
    txt = response.read()

    print 'downloading...' +  pdb_id + '_ba.pdb',

    f = tempfile.NamedTemporaryFile(delete=False)
    with open(f.name, 'w') as f:
        f.write(txt)

    with gzip.open(f.name, 'rb') as f:
        file_content = f.read()
        with open(pdb_id + '_ba.pdb', 'w') as ff:
            ff.write(file_content)
    print 'ok'


# main
if '__main__' == __name__:
    fn = 'input/image'
    print 'fn:', fn
    struc = StrucFile(fn)
    print ' pdb?:', struc.is_it_pdb()
    print ' # atoms:', struc.get_no_lines()

    fn = 'input/na.pdb'
    s = StrucFile(fn)
    print s.detect_molecule_type()
    #res = get_all_res(na)
    #print 'what is?', what_is(res)
    #print res
    print 'non standard:', s.check_res_if_std_na()
    print 'is protein:', s.detect_molecule_type()

    fn = 'input/prot.pdb'
    s = StrucFile(fn)
    print 'non standard:', s.check_res_if_std_prot()
    print 'is protein:',  s.detect_molecule_type()


    fn = 'input/rna-ru.pdb'
    s = StrucFile(fn)
    print 'non standard:', s.check_res_if_supid_rna()
    print 'is protein:', s.detect_molecule_type()

    fn = 'input/na_highAtomNum.pdb'
    print fn
    s = StrucFile(fn)
    s.renum_atoms()
    s.write('output/na_highAtomNum.pdb')

    fn = 'input/na_solvet_old_format.pdb'
    print fn
    s = StrucFile(fn)
    s.fix_op_atoms()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.write('output/na_solvet_old_format.pdb')

    fn = 'input/na_solvet_old_format.pdb'
    print fn
    s = StrucFile(fn)
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.write('output/na_solvet_old_format.pdb')

    #fn = 'input/na_solvet_old_format__.pdb'
    #s = StrucFile(fn)
    #s.fix_resn()
    #s.remove_hydrogen()
    #s.remove_ion()
    #s.remove_water()
    #s.renum_atoms()
    #s.fix_op_atoms()
    #s.write('output/na_solvet_old_format__.pdb')


    fn = 'input/1xjr.pdb'
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.renum_atoms()
    s.fix_op_atoms()
    s.write('output/1xjr.pdb')

    fn = 'input/decoy0165_amb.pdb'
    print fn
    s = StrucFile(fn)
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.renum_atoms()
    s.fix_O_in_UC()
    s.fix_op_atoms()
    s.write('output/decoy0165_amb_clx.pdb')

    fn = 'input/farna.pdb'
    print fn
    s = StrucFile(fn)
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.fix_op_atoms()
    s.renum_atoms()
    s.write('output/farna.pdb')

    fn = 'input/farna.pdb'
    print fn

    r = StrucFile(fn)
    print r.is_mol2()

    if True:
        print '================================================'
        print "input/1xjr_clx_fChimera_noIncludeNumbers.mol2"
        r = StrucFile("input/1xjr_clx_fChimera_noIncludeNumbers.mol2")
        print r.is_mol2()
        r.mol2toPDB('/tmp/x.pdb')

        r = StrucFile('/tmp/x.pdb')
        print r.get_report()
        r.fix_resn()
        r.remove_hydrogen()
        r.remove_ion()
        r.remove_water()
        r.fix_op_atoms()
        r.renum_atoms()
        r.fixU__to__U()
        r.write("output/1xjr_clx_fChimera_noIncludeNumbers.mol2")

    if True:
        r = StrucFile("input/2du3_prot_bound.mol2")
        print r.is_mol2()
        outfn = r.mol2toPDB()
        print r.get_report()

    print '================================================'
    fn = "input/3e5fA-nogtp_processed_zephyr.pdb"
    r = StrucFile(fn)
    print r.is_mol2()
    #outfn = r.mol2toPDB()
    print r.is_amber_like()
    print r.get_report()

    print r.get_preview()

    r.fix_resn()

    print r.get_preview()

    r.remove_hydrogen()
    r.remove_ion()
    r.remove_water()
    #renum_atoms(t, t)
    #fix_O_in_UC(t, t)
    #fix_op_atoms(t, t)
    r.write('output/3e5fA-nogtp_processed_zephyr.pdb')

    print
    fn = "input/1xjr_clx_charmm.pdb"
    print fn
    s = StrucFile(fn)
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.write('output/1xjr_clx_charmm.pdb')

    #renum_atoms(t, t)
    #fix_O_in_UC(t, t)
    #fix_op_atoms(t, t)

    print
    fn = "input/dna_fconvpdb_charmm22.pdb"
    print fn
    r = StrucFile(fn)
    r.get_preview()
    r.resn_as_dna()
    r.remove_hydrogen()
    r.remove_ion()
    r.remove_water()
    r.fix_resn()
    print r.get_head()
    print r.get_tail()
    print r.get_preview()
    r.write("output/dna_fconvpdb_charmm22.pdb")


    print
    fn = "input/1a9l_NMR_1_2_models.pdb"
    print fn
    r = StrucFile(fn)
    r.write("output/1a9l_NMR_1_2_models_lib.pdb")
    #r.get_text() # get #1 model
