#!/usr/bin/env python
#
# BaseRecognizer.py
#
# Class recognizing modified base residues.
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


from rna_tools.tools.mini_moderna3.moderna.analyze.MolecularGraph import AnnotatedMolecule
from rna_tools.tools.mini_moderna3.moderna.analyze.MolTopologies import read_nucleotide_topologies
from rna_tools.tools.mini_moderna3.moderna.util.Errors import BaseRecognitionError
from rna_tools.tools.mini_moderna3.moderna.Constants import AMINO, STANDARD_BASES, HETERO_GROUPS, PHOSPHATE_GROUP, \
    NUCLEOTIDE_ATOMS
from rna_tools.tools.mini_moderna3.moderna.Constants import MODIFICATION_TOPOLOGY_FILE
import re

MODIFIED_BASE_PATTERNS = read_nucleotide_topologies(MODIFICATION_TOPOLOGY_FILE)

# do not recognize amino acids
EXCLUDE_AMINO = True

DNA_BASES = ['dA', 'dC', 'dG', 'dT']

NAME_MATCHES = {
    'RMP': 'd5mpA', 
    'SMP': 'd5mpA', 
    'CMR':'d5mpC', 
    'A3A': 'alpha-dA', 
    'GAO': 'arabinoseG', 
    'CAR': 'arabinoseC', 
    'UAR': 'arabinoseU', 
    'A5O': 'arabinoseA', 
    'G25': 'GMP', 
    'A5L': 'dA', 
    'DA': 'dA', 
    
    }
    
    
class BaseRecognitionResult(object):
    
    def __init__(self, resi):
        self.resi = resi
        self.abbrev = ''
        self.mol = None
        self.subunits = []

    @property
    def amino_acid(self):
        """
        Returns amino acid three-letter code, if this is an amino acid,
        otherwise None.
        """
        short_name = self.resi.resname
        if short_name in AMINO:
            self.abbrev = short_name
            return short_name

    @property
    def standard_nucleotide(self):
        """
        Recognize nucleotides in a Bio.PDB.Residue instance that
        - have a name like A,G,C,T,U
        - only atoms with standard names in the PDB file
        - All atoms must be present, except the phosphate group.
        - The O2' atom is used to distinguish ribo- from desoxyribonucleotides.
        - all hydrogens are neglected
        returns one of A,G,C,T,U,dA,dG,dC,dT,dU, or False.
        """
        short_name = self.resi.resname.strip()
        atoms = set([])
        if short_name in STANDARD_BASES:
            desoxynucleotide = 'd'
            for atom in self.resi.child_list:
                atomname = atom.id.replace('*',"'")
                if atomname in PHOSPHATE_GROUP: continue
                elif atomname[0] == 'H': continue
                elif atomname == "O2'": desoxynucleotide = ''
                else:
                    atoms.add(atomname)
            if atoms==NUCLEOTIDE_ATOMS[short_name]:
                result = desoxynucleotide + short_name
                self.abbrev = result
                return result
            return False

    def has_forbidden_elements(self):
        FORBIDDEN = ['BR', 'CL', 'MN', 'MG',  'I', 'F', 'FE', 'CU', 'V', 'CR']
        for at in self.resi:
            name = re.sub("[\d\s']", '', at.fullname)
            name = name.upper()
            if name in FORBIDDEN:
                return True
        
    def check_forbidden_elements(self):
        """catch heavy atom groups."""
        if self.has_forbidden_elements():
            if self.abbrev in STANDARD_BASES + DNA_BASES:
                self.abbrev = '?' + self.abbrev
            elif self.abbrev == '':
                raise BaseRecognitionError('Strange element occured in residue recognized as : %s'%self.abbrev)

    def identify_phosphates(self, tags, resn):
        """Distinguishes ATP, ADP, AMP, GTP, ... and other phosphate extensions from standard nucleotides."""
        restype = NAME_MATCHES.get(resn)
        if restype in tags:
            return restype

        for t in tags:
            if t.endswith('_mod_phos'):
                return t
            elif resn == t or 'd'+resn ==t:
                restype = t
            elif not restype and t in STANDARD_BASES + DNA_BASES:
                restype = t
            
        #if not restype and len(tags)>0 and tags[0] in STANDARD_BASES + DNA_BASES:
        #    restype = tags[0]
        return restype

    def run_topology_matching(self):
        """Running the subgraph matching algorithm."""
        self.mol = AnnotatedMolecule()
        self.mol.parse_resi(self.resi)
        self.subunits = set([su for su, atom in self.mol.detect_patterns(MODIFIED_BASE_PATTERNS)])
        #print self.subunits

    def identify_molecule(self):
        """returns an abbreviated name for the molecule"""
        restype = ""
        resn = re.sub('\s','',self.resi.resname.strip())
        sub2 = [x for x in self.subunits if x not in ('phosphate','desoxyriboside','riboside','purine','pyrimidine')]
        sub2.sort()

        restype = self.identify_phosphates(sub2, resn)
        # decide about some difficult cases
        if not restype:
            if len(sub2) == 1:
                if sub2[0] in STANDARD_BASES: # the easy ones
                    restype = sub2[0]
                elif sub2[0]=='phosphate': restype = '' # phosphate only
                else: restype = '['+sub2[0]+']' # put modified bases in [x]
            elif sub2 == ['m5D', 'm5U']:
                if re.search('D',resn): restype = '[m5D]'
                else: restype = prefix+'[m5U]'
            elif sub2 == ['D', 'U']:
                if re.search('D',resn): restype = '[D]'
                else: restype = prefix+'U'
            elif sub2 == ['T','m5D','m5U'] or sub2 == ['T','m5U']:
                if resn == 'T': restype = 'T'
                elif resn == 'U': restype = 'U'
                else: restype = '[m5U]'
            elif sub2 == ['galQtRNA', 'manQtRNA']:
                if re.search('g',resn,re.IGNORECASE):
                    restype = '[galQtRNA]'
                else:
                    restype ='[manQtRNA]'
            else:
                restype = '<%s:%s>'%(resn,str(list(self.subunits)))
        
            if restype[0]=='<':
                if restype[1:-1] in HETERO_GROUPS:
                    restype = ''
                else:
                    # probably some other hetero group,
                    raise BaseRecognitionError('Unknown Residue:%s'%restype)
        
        restype = re.sub('[\[\]]','',restype)
        self.abbrev = restype
        return restype



class BaseRecognizer(object):
    """
    Assigns a name to a residue - if necessary by analyzing
    the detailed topolgy of atoms.
    """
    def identify_resi(self,resi):
        """
        Recognizes what kind of residue there is and returns an abbreviation.
        Takes a Bio.PDB.Residue instance.
        
        First, assignment by atom and residue names will be done.
        - Names must be A,G,C,T,U

        Second, the structure will be converted to .mol format and
        examined by the pattern matching procedure.
        """
        result = BaseRecognitionResult(resi)
        if result.amino_acid or result.standard_nucleotide:
            return result.abbrev

        # Now, this residue might be a modified base,
        result.run_topology_matching()
        restype = result.identify_molecule()
        if not result.abbrev:
            raise BaseRecognitionError('Unknown Residue:%s'%restype)
            
        result.check_forbidden_elements()
        return result.abbrev
        
