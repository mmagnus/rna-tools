#!/usr/bin/python

from SingleLineUtils import get_res_code, get_res_num, get_atom_code, \
    set_atom_code
import PDBFile
import os

if __name__ == '__main__':
    fn = '1xjrA_M1.pdb'
    pdb_file = PDBFile.PDBFile(pdb_path='test' + os.sep + fn)
    print((pdb_file.pdb_lines[0]))
    pdb_file.pedantic_pdb()
    print((pdb_file.fixes))
    print((type(pdb_file.pdb_lines)))
    print(('14:', pdb_file.pdb_lines[0]))
    pdb_file.remove_non_atoms()
    print((pdb_file.pdb_lines[0]))
    pdb_file.check_and_add_P_at_start()
    pdb_file.check_and_add_P_at_start()
    print((pdb_file.fixes))
    pdb_file.save('tmp.pdb')
