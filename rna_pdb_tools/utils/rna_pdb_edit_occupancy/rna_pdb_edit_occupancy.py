from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.Atom import PDBConstructionWarning

import sys
import warnings
warnings.simplefilter('ignore', PDBConstructionWarning)

import re
import string

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

if __name__ == '__main__':
    #freeze = "A:1-10,15,25-30 ; B:1-10" # error
    #freeze = "A:1-10,15,25-300" # error
    #freeze = "A:1-2,15,25-24" # error
    freeze = "A:1-2, 4" # error
    pdb = "test_data/3w3s_homologymodel.pdb"
    pdb_out = 'out.pdb'
    edit_occupancy_of_pdb(freeze, pdb, pdb_out)


    
