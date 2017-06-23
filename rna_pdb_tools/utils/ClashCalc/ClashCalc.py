#!/usr/bin/python

from Bio.PDB import NeighborSearch, PDBParser, Selection, Atom
from numpy import array

def check_clash(str_name, v=True):
        """check_clash, fract of clashes!

        if zero contacts then error -> fix ->

        Problem, contacts, str_name: 311 505 na-prot_13536.pdb
        Sterical clashes  0.615841584158

        c is counter
        """
        print(str_name)
        structure = open(str_name)
        #model = structure[0]
        atoms_A = []
        atoms_B = []
        for line in structure.readlines():
            if line[:4] == "ATOM":
                #print line
                at_nam = line[12:16].strip()
                coor = [float(line[30:38]),float(line[38:46]), float(line[46:54])]	
                at = Atom.Atom(at_nam,coor,0.0,1.0,' ',at_nam,1,at_nam[0])
                if line[21] == "A":
                    atoms_A.append(at)
                elif line[21] == "B":
                    atoms_B.append(at)
                else: pass
        #atoms_B = Selection.unfold_entities(structure[0]['B'], 'A')
        #print len(atoms_A), len(atoms_B)
        if len(atoms_A) > len(atoms_B):
            less = atoms_B
            more = atoms_A
        else: 
            less = atoms_A
            more = atoms_B
        problem = 0
        contacts = 0 
        ns=NeighborSearch(more)
        for at in less:
             neighbors=ns.search(array(at.get_coord()),2.0,'A')
             if neighbors != []:
                 problem +=1
                 contacts +=1
             else:
                 neighbors1=ns.search(array(at.get_coord()),4.0,'A')
                 if neighbors1 != []:
                     contacts +=1
        if v:
                print('problem:', float(problem))
                print('contacts:', float(contacts))
        try:
            fract = float(problem)/float(contacts)
        except ZeroDivisionError:
            fract = problem # or skip this structure
            print('ZeroDivison -- skip:', problem, contacts, str_name)
            return fract

        #print 'Contacts, str_name:', problem, contacts, str_name, "Sterical clashes ", fract
        return fract
    
if __name__ == '__main__':
    print(check_clash('test_data/no_clash.pdb'))
    print(check_clash('test_data/super_clash.pdb'))
    print(check_clash('test_data/prot-na_2392.pdb'))
