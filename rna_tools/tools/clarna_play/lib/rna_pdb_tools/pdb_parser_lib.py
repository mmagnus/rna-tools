#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from __future__ import print_function
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

ignore_op3 = False

def get_version(currfn='', verbose=False): #dupa
    """Get version of the tool based on state of the git repository.
    Return version. 
    If currfn is empty, then the path is '.'. Hmm.. I think it will work. We will see.
    The version is not printed!
    https://github.com/m4rx9/curr_version/"""
    from subprocess import getoutput

    if currfn == '':
        path = '.'
    else:
        path = os.path.dirname(currfn)
    if verbose: print('get_version::path', path)
    if os.path.islink(currfn):#path + os.sep + os.path.basename(__file__)):
        path = os.path.dirname(os.readlink(path + os.sep + os.path.basename(currfn)))
    if not path: path = '.'
    if verbose: print('get_version::path2', path)
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
            print('pybel is needed for mol2 to pdb convertion')
            #sys.exit(1)
            sys.exit(0)

        if not outfn:
            outfn = self.fn.replace('.mol2', '.pdb')
            
        for mol in pybel.readfile("mol2", self.fn):
            mol.write("pdb", outfn, overwrite=True)

        print('outfn: ', outfn)
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
        """
        You get `chains` such as:
        OrderedDict([('A', {'header': 'A:1-47', 'seq': 'CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAU'}), ('B', {'header': 'B:48-58', 'seq': 'AUCAGGUGCAA'})])"""
        
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
                    if chain_curr in chains:                
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
        for c in list(chains.keys()):
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
                    if chain_curr in chains:                
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
        for c in list(chains.keys()):
            seq += chains[c]['header'].replace(':', ' ').replace('-', ' ') + ' '
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
        print('fixU__to__U OK')
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

        print('resn_as_dna')
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
        lines = []
        for l in self.lines:
            nl = l.replace('*', '\'')
            nl = nl.replace('O1P', 'OP1')
            nl = nl.replace('O2P', 'OP2')
            nl = nl.replace('O3P', 'OP3')
            lines.append(nl)
        self.lines = lines

    def get_report(self):
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
            print('Write %s' % outfn)

    def get_rnapuzzle_ready(self, renumber_residues=True):
        """Get rnapuzzle ready structure.
        Submission format @http://ahsoka.u-strasbg.fr/rnapuzzles/

        Does:
        - keep only given atoms,
        - renumber residues from 1, if renumber_residues=True (by default)
        """
        try:
            from Bio import PDB
            from Bio.PDB import PDBIO
        except:
            sys.exit('Error: Install biopython to use this function (pip biopython)')

        import copy

        G_ATOMS = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'', 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4']
        A_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split()
        U_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split()
        C_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split()

        ftmp = '/tmp/out.pdb'
        self.write(ftmp,v=False)

        parser = PDB.PDBParser()
        struct = parser.get_structure('', ftmp)
        model = struct[0]

        s2 = PDB.Structure.Structure(struct.id)
        m2 = PDB.Model.Model(model.id)

        chains2 = []

        missing = []
        for chain in model.get_list():
            res = [] 
            for r in chain:
                res.append(r)

            res = copy.copy(res)

            c2 = PDB.Chain.Chain(chain.id)        

            c = 1  # new chain, goes from 1 !!!
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
                if str(r.get_resname()).strip() == "G":

                    for an in G_ATOMS:
                        try:
                            r2.add(r[an])
                        except KeyError:
                            #print 'Missing:', an, r, ' new resi', c
                            missing.append([an, chain.id, r, c])
                    c2.add(r2)

                elif str(r.get_resname()).strip() == "A":
                    for an in A_ATOMS:
                        try:
                            r2.add(r[an])
                        except KeyError:
                            #print 'Missing:', an, r, ' new resi', c
                            missing.append([an, chain.id, r, c])
                    c2.add(r2)

                elif str(r.get_resname()).strip() == "C":
                    for an in C_ATOMS:
                        try:
                            r2.add(r[an])
                        except:
                            #print 'Missing:', an, r, ' new resi', c
                            missing.append([an, chain.id, r, c])
                    c2.add(r2)

                elif str(r.get_resname()).strip() == "U":
                    for an in U_ATOMS:
                        try:
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
        #fout = fn.replace('.pdb', '_fx.pdb')
        fout = '/tmp/outout.pdb' # hack
        io.save(fout)
        
        if missing:
            print('REMARK 000 Missing atoms:')
            for i in missing:
                print('REMARK 000  +', i[0], i[1], i[2], 'residue #', i[3])
            #raise Exception('Missing atoms in %s' % self.fn)
        s = StrucFile(fout)
        self.lines = s.lines


    def get_simrna_ready(self,  renumber_residues=True):
        """Get simrna_ready .. 

        - take only first model,
        - renumber residues if renumber_residues=True

        .. warning:: requires: Biopython"""
        try:
            from Bio import PDB
            from Bio.PDB import PDBIO
        except:
            sys.exit('Error: Install biopython to use this function (pip biopython)')

        import warnings
        
        warnings.filterwarnings('ignore', '.*Invalid or missing.*',)
        warnings.filterwarnings('ignore', '.*with given element *',)
        
        import copy

        G_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split()
        A_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split()
        U_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split()
        C_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split()

        ftmp = '/tmp/out.pdb'
        self.write(ftmp,v=False)

        parser = PDB.PDBParser()
        struct = parser.get_structure('', ftmp)
        model = struct[0]

        s2 = PDB.Structure.Structure(struct.id)
        m2 = PDB.Model.Model(model.id)

        chains2 = []

        missing = []
        
        for chain in model.get_list():
            res = [] 
            for r in chain:
                res.append(r)

            res = copy.copy(res)

            c2 = PDB.Chain.Chain(chain.id)        

            c = 1  # new chain, goes from 1 if renumber True
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
                if c == 1:
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

                    if p_missing:
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
        #fout = fn.replace('.pdb', '_fx.pdb')
        fout = '/tmp/outout.pdb' # hack
        io.save(fout)
        
        if missing:
            print('REMARK 000 Missing atoms:')
            for i in missing:
                print('REMARK 000  +', i[0], i[1], i[2], 'residue #', i[3])
            #raise Exception('Missing atoms in %s' % self.fn)
        s = StrucFile(fout)
        self.lines = s.lines

    def edit_occupancy_of_pdb(txt, pdb, pdb_out,v=False):
        """Make all atoms 1 (flexi) and then set occupancy 0 for seletected atoms.
        Return False if error. True if OK
        """
        struc = PDB.PDBParser().get_structure('struc', pdb)

        txt = txt.replace(' ','')
        if v:print(txt)
        l = re.split('[,:;]', txt)
        if v:print(l) 

        for s in struc:
            for c in s:
                for r in c:
                    for a in r:
                        a.set_occupancy(1)  # make it flaxi

        for i in l: # ['A', '1-10', '15', '25-30', 'B', '1-10']

            if i in string.ascii_letters:
                if v:print('chain', i)
                chain_curr = i
                continue

            if i.find('-') > -1:
                start, ends = i.split('-')
                if start > ends:
                    print('Error: range start > end ' + i, file=sys.stderr)
                    return False
                index = list(range(int(start), int(ends)+1))
            else:
                index=[int(i)]

            for i in index:
                # change b_factor
                try:
                    atoms = struc[0][chain_curr][i]
                except KeyError:
                    if i == chain_curr:
                        print('Error: Chain ' + chain_curr + ' not found in the PDB structure', file=sys.stderr)
                    else:
                        print('Error: Residue ' + chain_curr + ':' + str(i) + ' found in the PDB structure', file=sys.stderr)
                        return False
                for a in atoms:
                    a.set_occupancy(0)

        io = PDBIO()
        io.set_structure(struc)
        io.save(pdb_out)
        print('Saved ', pdb_out)
        return True

# main
if '__main__' == __name__:
    fn = 'input/image'
    print('fn:', fn)
    struc = StrucFile(fn)
    print(' pdb?:', struc.is_it_pdb())
    print(' # atoms:', struc.get_no_lines())

    fn = 'input/na.pdb'
    s = StrucFile(fn)
    print(s.detect_molecule_type())
    #res = get_all_res(na)
    #print 'what is?', what_is(res)
    #print res
    print('non standard:', s.check_res_if_std_na())
    print('is protein:', s.detect_molecule_type())

    fn = 'input/prot.pdb'
    s = StrucFile(fn)
    print('non standard:', s.check_res_if_std_prot())
    print('is protein:',  s.detect_molecule_type())


    fn = 'input/rna-ru.pdb'
    s = StrucFile(fn)
    print('non standard:', s.check_res_if_supid_rna())
    print('is protein:', s.detect_molecule_type())

    fn = 'input/na_highAtomNum.pdb'
    print(fn)
    s = StrucFile(fn)
    s.renum_atoms()
    s.write('output/na_highAtomNum.pdb')

    fn = 'input/na_solvet_old_format.pdb'
    print(fn)
    s = StrucFile(fn)
    s.fix_op_atoms()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.write('output/na_solvet_old_format.pdb')

    fn = 'input/na_solvet_old_format.pdb'
    print(fn)
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
    print(fn)
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
    print(fn)
    s = StrucFile(fn)
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.fix_op_atoms()
    s.renum_atoms()
    s.write('output/farna.pdb')

    fn = 'input/farna.pdb'
    print(fn)

    r = StrucFile(fn)
    print(r.is_mol2())

    if True:
        print('================================================')
        print("input/1xjr_clx_fChimera_noIncludeNumbers.mol2")
        r = StrucFile("input/1xjr_clx_fChimera_noIncludeNumbers.mol2")
        print(r.is_mol2())
        r.mol2toPDB('/tmp/x.pdb')

        r = StrucFile('/tmp/x.pdb')
        print(r.get_report())
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
        print(r.is_mol2())
        outfn = r.mol2toPDB()
        print(r.get_report())

    print('================================================')
    fn = "input/3e5fA-nogtp_processed_zephyr.pdb"
    r = StrucFile(fn)
    print(r.is_mol2())
    #outfn = r.mol2toPDB()
    print(r.is_amber_like())
    print(r.get_report())

    print(r.get_preview())

    r.fix_resn()

    print(r.get_preview())

    r.remove_hydrogen()
    r.remove_ion()
    r.remove_water()
    #renum_atoms(t, t)
    #fix_O_in_UC(t, t)
    #fix_op_atoms(t, t)
    r.write('output/3e5fA-nogtp_processed_zephyr.pdb')

    print()
    fn = "input/1xjr_clx_charmm.pdb"
    print(fn)
    s = StrucFile(fn)
    s.fix_resn()
    s.remove_hydrogen()
    s.remove_ion()
    s.remove_water()
    s.write('output/1xjr_clx_charmm.pdb')

    #renum_atoms(t, t)
    #fix_O_in_UC(t, t)
    #fix_op_atoms(t, t)

    print()
    fn = "input/dna_fconvpdb_charmm22.pdb"
    print(fn)
    r = StrucFile(fn)
    r.get_preview()
    r.resn_as_dna()
    r.remove_hydrogen()
    r.remove_ion()
    r.remove_water()
    r.fix_resn()
    print(r.get_head())
    print(r.get_tail())
    print(r.get_preview())
    r.write("output/dna_fconvpdb_charmm22.pdb")


    print()
    fn = "input/1a9l_NMR_1_2_models.pdb"
    print(fn)
    r = StrucFile(fn)
    r.write("output/1a9l_NMR_1_2_models_lib.pdb")
    #r.get_text() # get #1 model
