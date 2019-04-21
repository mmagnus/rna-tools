#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rna_tools.tools.pdb_formatix.SingleLineUtils import get_res_code, get_res_num, get_atom_code, \
    set_atom_code, set_line_bfactor
import re


class PDBFile(object):
    """Class for holding data from a PDB file and modifying it.
    """

    # find 'ATOM' lines in a PDB file
    ATOM_LINE_PATTERN = re.compile('^ATOM')

    # dictionary for changing residue names from 3 letters to 1 letter
    RES3TO1 = {
        'GUA': 'G',
        'URI': 'U',
        'CYT': 'C',
        'ADE': 'A',

        'RG': 'G',
        'RU': 'U',
        'RC': 'C',
        'RA': 'A',

        'R3G': 'G',
        'R3U': 'U',
        'R3C': 'C',
        'R3A': 'A',

        'R5G': 'G',
        'R5U': 'U',
        'R5C': 'C',
        'R5A': 'A',

        'RG3': 'G',
        'RU3': 'U',
        'RC3': 'C',
        'RA3': 'A',

        'RG5': 'G',
        'RU5': 'U',
        'RC5': 'C',
        'RA5': 'A',

        'RGU': 'G',
        'URA': 'U',
        'RCY': 'C',
        'RAD': 'A',

        # just in case
        'G': 'G',
        'U': 'U',
        'C': 'C',
        'A': 'A',
    }

    AMINOACID_CODES = ("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly",
                       "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr",
                       "Trp", "Tyr", "Val",
                       "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                       "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
                       "TRP", "TYR", "VAL"
                       )

    def __init__(self, pdb_string=None, pdb_path=None, pdb_handle=None, verbose=False):
        """Constructor, should get exactly one of the arguments:
             * pdb_string = PDB as a string
             * pdb_path = string with path to a PDB file
             * pdb_handle = handle to a PDB file
        """
        self.verbose = verbose

        if len([x for x in [pdb_string, pdb_path, pdb_handle] if x is not None]) > 1:
            print('You should provide at most one source for PDB file')
            raise Exception
        input_string = ''
        if pdb_string is not None:
            input_string = pdb_string
        elif pdb_path is not None:
            with open(pdb_path) as f:
                input_string = f.read()
        elif pdb_handle is not None:
            input_string = pdb_handle.read()
        self.pdb_string = input_string
        self.pdb_lines = self.pdb_string.split('\n')
        self.fixes = []

    def save(self, file_path):
        """Save current PDB to disk

        Arguments:
          * file_path = path where it will be saved
        """
        with open(file_path, 'w') as f:
            f.write(self.pdb_string)

    def _apply_fix(self, fix_name, result_string, result_lines=None):
        """Helper function for applying fixes and saving information about them.

        Arguments:
          * fix_name = string that will be added to self.pdb_fixes, unless
            result_string is the same as pdb_string
          * result_string = string after applying this fix
          * result_lines = optional, list of lines, use this argument if you
            already have such list and don't want to lose time on splitting
            result_string
        """
        if self.pdb_string != result_string:
            self.fixes.append(fix_name)
            self.pdb_string = result_string
            if result_lines is not None:
                self.pdb_lines = result_lines
            else:
                self.pdb_lines = self.pdb_string.split('\n')

    def set_string(self, pdb_string):
        """Change PDB string stored in this instance of the class

        Arguments:
          * pdb_string = new PDB string
        """
        self.pdb_string = pdb_string
        self.pdb_lines = pdb_string.split('\n')

    def _get_atom_lines(self):
        """Get only lines with ATOM information
        """
        return [l for l in self.pdb_lines if re.match(self.ATOM_LINE_PATTERN, l)]

    def validate_pdb(self):
        """Check if file is a PDB structure

        Output:
          * True if it is a PDB, False otherwise
        """
        atom_lines = self._get_atom_lines()
        if not len(atom_lines):
            return False
        return True

    def detect_proteins(self):
        """Check if there are any aminoacid fragments in the file

        Output:
          * True if there are some, False otherwise
        """
        for l in self.pdb_lines:
            if get_res_code(l) in self.AMINOACID_CODES:
                return True
        return False

    def seq_from_pdb(self):
        """Extract sequence from a PDB and return it as a string.

        Output:
          * sequence, returned as a string
        """
        atoms = [l.split() for l in self._get_atom_lines()]
        if atoms[0][3][0] != 'r':  # it is 'r' if it's a ROSETTA PDB
            seq = [self.RES3TO1[atoms[0][3]]]
        else:
            seq = [atoms[0][3][1]]
        for a in range(1, len(atoms)):
            # atom number is different than previous one
            if atoms[a][5] != atoms[a - 1][5]:
                if atoms[a][3][0] != 'r':  # check for ROSETTA PDB
                    seq.append(self.RES3TO1[atoms[a][3]])
                else:
                    seq.append(atoms[a][3][1])
        return ''.join(seq)

    def seq_from_amber_like_pdb(self):
        """Extract sequence from a PDB and return it as a string - use it for amber-like files.

        Output:
          * sequence, returned as a string, such as: RG5 RC RU RG RG RG RC RG RC RA RG RG3 RC5 RC RU RG RA RC RG RG RU RA RC RA RG RC3
        """
        atoms = [l.split() for l in self._get_atom_lines()]
        seq = []
        seq.append(atoms[0][3])
        for a in range(1, len(atoms)):
            # atom number is different than previous one
            if self.verbose:
                print((atoms[a][5], atoms[a - 1][5]))
            if atoms[a][5] != atoms[a - 1][5]:
                seq.append(atoms[a][3])
        return ' '.join(seq)

    def get_fasta(self, name='seq', lowercase=False):
        """Format sequence in FASTA format, with a header and lines split
        at 80 characters.

        Arguments:
          * name = name of the sequence (it's put in the header line)

        Output:
          * FASTA returned as a string
        """
        seq = self.seq_from_pdb()
        if lowercase:
            seq = seq.lower()
        result = ['>' + name]
        result.extend([seq[i:i + 80] for i in range(0, len(seq), 80)])
        result.append('')
        return '\n'.join(result)

    def remove_non_atoms(self):
        """Remove all lines that are not ATOMs
        """
        result = self._get_atom_lines()
        self._apply_fix('Removed non-atom lines', '\n'.join(result), result)

    def _check_resname_3(self):
        """Check if PDB uses 3 or 1 letter residue names.

        Output:
          * bool, True if 3 letters, False otherwise
        """
        for l in self.pdb_lines:
            if l.startswith('ATOM'):
                return 3 == len(l.split()[3])

    def _resname_3to1(self):
        """Convert residue names in PDB file from 3 letters to 1 letter.
        """
        result = []
        for l in self.pdb_lines:
            if l.startswith('ATOM'):
                long_name = l.split()[3]
                try:
                    short_name = self.RES3TO1[long_name]
                except KeyError:
                    short_name = 'X'
                result.append(l[:17] + '  ' + short_name + l[20:])
            else:
                result.append(l)
        self._apply_fix('resname_3to1', '\n'.join(result), result)

    def resname_check_and_3to1(self):
        """Check if resnames are 3 letter long and if so convert them to 1 letter.
        """
        if self._check_resname_3():
            try:
                self._resname_3to1()
            except:
                print('Conversion to 1-letter residue names failed')

    def terminate_chains(self):
        """Add 'TER' at the end of chain if none 'TER's are found in PDB.
        """
        for l in self.pdb_lines:
            if l.startswith('TER'):
                return
        else:
            result = []
            chain_started = False
            ter_added = False
            for l in self.pdb_lines:
                if 'ATOM' in l:
                    chain_started = True
                else:
                    if chain_started:
                        result.append('TER')
                        ter_added = True
                        chain_started = False
                result.append(l)
        if not ter_added:
            result.append('TER')
        self._apply_fix('terminate_chains', '\n'.join(result), result)

    def _split_by_ters(self, separator='TER\n'):
        """Split a PDB string by TER lines or other separator

        Arguments:
          * separator = optional, separator dividing structures, default is 'TER\n'
        """
        return [i + separator for i in self.pdb_string.split(separator)]

    def remove_short_chains(self, threshold=9):
        """Remove chains that have chains with a small number of atoms.

        Arguments:
          * threshold = if chain has more atoms, it stays in the result
        """
        chains = self.pdb_string.split('TER\n')
        good_chains = [c for c in chains if len(c.split('\nATOM')) > threshold]
        self._apply_fix('remove_short_chains', '\n'.join(good_chains))

    def check_and_add_P_at_start(self):
        residues = {}
        for l in self.pdb_lines:
            res_num = get_res_num(l)
            if res_num in residues:
                residues[res_num].append(l)
            else:
                residues[res_num] = [l]
        first_residue = residues[min(residues.keys())]
        # check P
        has_p = False
        for l in first_residue:
            if get_atom_code(l) == 'P':
                has_p = True
        if not has_p:
            # add P
            corrected_residue = []
            for l in first_residue:
                if get_atom_code(l) == 'O5\'' or get_atom_code(l) == 'O5*':
                    corrected_residue.append(set_atom_code(l, 'P'))
                else:
                    corrected_residue.append(l)
            residues[min(residues.keys())] = corrected_residue
            residues = ['\n'.join(residues[r]) for r in residues]
            self._apply_fix('add_P_at_start', '\n'.join(residues))

    def set_residues_bfactor(self, bfactors):
        """Set B-factor to any value you want

        Arguments:
          * bfactors = list of B-factors that will be set, should be of same
            length as nucleotide sequence
        """
        ans = []
        for l in self.pdb_lines:
            if l.startswith('ATOM'):
                res_num = get_res_num(l)
                try:
                    ans.append(set_line_bfactor(l, bfactors[res_num - 1]))
                except IndexError:
                    pass
            else:
                ans.append(l)
        self.pdb_lines = ans
        self.pdb_string = '\n'.join(ans)

    def pedantic_pdb(self):
        """Do everything that's possible to fix PDB: 3-to-1, no HETATMs, TERs etc.
        """
        self.resname_check_and_3to1()
        self.terminate_chains()
        self.remove_short_chains()

    def count_models(self):
        """Count models in the PDB file

        Output:
          * number of models as an int
        """
        model_num = 0
        for l in self.pdb_lines:
            if l.startswith('MODEL'):
                model_num += 1
        return model_num

    def get_model(self, model_num):
        """Get n-th model from the file

        Arguments:
          * model_num = number of model to get, starts with 0
        """
        model_borders = list(zip(
            [i[0] for i in enumerate(self.pdb_lines) if i[1].startswith('MODEL')],
            [i[0] for i in enumerate(self.pdb_lines) if i[1].startswith('ENDMDL')]
        ))
        result = self.pdb_lines[:model_borders[0][0]]
        result.extend(self.pdb_lines[model_borders[model_num][0]:model_borders[model_num][1] + 1])
        result.extend(self.pdb_lines[model_borders[-1][1] + 1:])
        self._apply_fix('get_model', '\n'.join(result), result)

    def check_and_get_first_model(self):
        """Check if there are more than one models and get only the first one
        """
        if self.count_models() >= 2:
            self.get_model(0)
