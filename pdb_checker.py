#!/usr

#from formatlib.PDBFile import PDBFile

AMINOACID_CODES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
            "TRP", "TYR", "VAL"]
RES = ['DA', 'DG', 'DT', 'DC']
RES += ['A', 'G', 'U', 'C']

def check_no_lines(pl):
    """Check number of lines of pdb file"""
    c = 0
    for l in pl:
        if l.startswith('ATOM'):
            c = c +1
    return c

def get_all_res(pl):
    res = []
    for l in pl:
        if l.startswith('ATOM'):
            r = l[17:20].strip().upper()
            if r not in res:
                res.append(r)
    return res

def check_res_if_std_na(res):
    wrong = []
    for r in res:
        if r not in RES:
            wrong.append(r)
    return wrong

def check_res_if_std_prot(res):
    wrong = []
    for r in res:
        if r not in AMINOACID_CODES:
            wrong.append(r)
    return wrong

def has_atom_line(na):
    for i in na:
        if i.startswith('ATOM'):
            return True
    return False

def get_atom_lines(fn):
    lines = open(fn).read().split('\n')
    al = []
    for l in lines:
        if l.startswith('ATOM'):
            al.append(l)
    return al


def check_res_if_std_dna(res):
    wrong_res = []
    for r in res:
        if r.upper().strip() in ['A', 'T', 'C', 'G']:
            if r not in wrong_res:
                wrong_res.append(r)
    return wrong_res

def check_res_if_supid_rna(res):
    wrong_res = []
    for r in res:
        if r.upper().strip() in ['RC', 'RU', 'RA', 'RG', 'RT']:
            if r not in wrong_res:
                wrong_res.append(r)
    return wrong_res

def is_rna(res):
    for r in res:
        
        if r.upper().strip() in ['RC', 'RU', 'RA', 'RG', 'RT']:
            if r not in wrong_res:
                wrong_res.append(r)
    return wrong_res

def renum_atoms(fn,out):
    c = 1
    ntxt = ''
    for l in open(fn).read().split('\n'):
        if l.startswith('ATOM') or l.startswith('HETATM') :
            print l
            nl = l[:6] + str(c).rjust(5) + l[11:]
            print nl
            c += 1
            ntxt += nl + '\n'

    o = open(out, 'w')
    o.write(ntxt)
    o.close()

if '__main__' == __name__:
    fn = 'test_data/image'
    print 'fn:', fn
    na = get_atom_lines(fn)
    print 'has_atom_line:', has_atom_line(na)

    fn = 'test_data/na.pdb'
    na = get_atom_lines(fn)
    print has_atom_line(na)

    c = check_no_lines(na)
    print c
    res = get_all_res(na)
    #print res
    print 'non standard:', check_res_if_std_na(res)
    #print 'is protein:', is_protein(res)

    fn = 'test_data/prot.pdb'
    prot = open(fn).read().split('\n')
    res = get_all_res(prot)
    print 'non standard:', check_res_if_std_prot(res)   
    #print 'is protein:', is_protein(res)


    fn = 'test_data/rna-ru.pdb'
    prot = open(fn).read().split('\n')
    res = get_all_res(prot)
    print 'non standard:', check_res_if_supid_rna(res)
    #print 'is protein:', is_protein(res)

    fn = 'test_data/na_highAtomNum.pdb'
    renum_atoms(fn, '/tmp/out.pdb')
