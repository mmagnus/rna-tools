#!/usr/bin/env python

#from formatlib.PDBFile import PDBFile

AMINOACID_CODES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
            "TRP", "TYR", "VAL"]
RES = ['DA', 'DG', 'DT', 'DC']
RES += ['A', 'G', 'U', 'C']

RESS = ['A', 'C', 'G', 'U', 'ADE', 'CYT', 'GUA', 'URY', 'URI', 'U34', 'U31', 'C31', '4SU', 'H2U', 'QUO', 'G7M', '5MU', '5MC', 'PSU', '2MG', '1MG', '1MA', 'M2G', '5BU', 'FHU', 'FMU', 'IU', 'OMG', 'OMC', 'OMU', 'A2M', 'A23', 'CCC', 'I'] + ['RC', 'RU', 'RA', 'RG', 'RT']
#DNA = ['DA', 'DG', 'DT', 'DC']
#RNA = ['A', 'G', 'U', 'C']

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

def what_is(res):
    """give it a score"""
    aa = []
    na = []

    for r in res:
        if r in AMINOACID_CODES:
            aa.append(r)
        if r in RESS:
            na.append(r)            

    aa = float(len(aa)) / len(res)
    na = float(len(na)) / len(res)

    if aa == 0 and na == 0:
        return 'error'
    if aa > na:
        return '>protein< vs na', aa, na
    else:
        return 'protein vs >na<', aa, na


def is_rna(res):
    wrong = []
    for r in res:
        if r.upper().strip() in ['RC', 'RU', 'RA', 'RG', 'RT']:
            if r not in wrong_res:
                wrong_res.append(r)
    return wrong_res


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
    """Renum atoms, from 1 to X for line; ATOM/HETATM"""
    c = 1
    ntxt = ''
    i = open(fn)
    txt = i.read().split('\n')
    for l in txt:
        if l.startswith('ATOM') or l.startswith('HETATM') :
            #print l
            nl = l[:6] + str(c).rjust(5) + l[11:]
            #print nl
            c += 1
            ntxt += nl + '\n'

        if l.startswith("END") or l.startswith("TER"):
            ntxt += l + '\n'

    i.close()
    o = open(out, 'w')
    o.write(ntxt)
    o.close()

def fix_op_atoms(fn,out):
    c = 1
    ntxt = ''
    i = open(fn)
    txt = i.read().split('\n')
    for l in txt:
        if l.startswith('ATOM') or l.startswith('HETATM') :
            l = l.replace('*', '\'')
            nl = l.replace('O1P', 'OP1')
            nl = nl.replace('O2P', 'OP2')
            nl = nl.replace('O3P', 'OP3')

            ntxt += nl + '\n'

        if l.startswith("END") or l.startswith("TER"):
            ntxt += l + '\n'

    i.close()
    o = open(out, 'w')
    o.write(ntxt)
    o.close()



def fix_O_in_UC(fn, out):
    """.. warning: remove RU names before using this function"""

    pdb_fn = fn
    f = open(pdb_fn)

    ntxt = ''

    for l in f:
        #if l[12:16].strip() in 
        #if l[12:16].strip().startswith('H'):
        nl = l.replace('O     U',
                       'O2    U')
        nl =nl.replace('O     C',
                       'O2    C')
        ntxt += nl
    f.close()

    f = open(out, 'w')
    f.write(ntxt)
    f.close()


def remove_hydrogen(fn, out):
    pdb_fn = fn
    f = open(pdb_fn)

    ntxt = ''
    hydrogen_names = ["H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", "H3", "H5", "H6", "H5T", "H41", "1H5'", 
                      "2H5'", "HO2'", "1H4", "2H4", "1H2", "2H2", "H1", "H8", "H2", "1H6", "2H6",
                      "HO5'", "H21", "H22", "H61", "H62", "H42", "HO3'", "1H2'", "2HO'", "HO'2", "H2'1" , "HO'2", "HO'2",
                      "H2", "H2'1", "H1", "H2", "1H5*","2H5*", "H4*", "H3*", "H1*", "1H2*", "2HO*", "1H2", "2H2", "1H4", "2H4", "1H6", "2H6", "H1", "H2", "H3", "H5", "H6", "H8", "H5'1", "H5'2"]

    for l in f:
        if l[77:79].strip() == 'H':
            continue
        if l[12:16].strip() in hydrogen_names:
        #if l[12:16].strip().startswith('H'):
            continue
        else:
            #print l[12:16]
            ntxt += l
    f.close()

    f = open(out, 'w')
    f.write(ntxt)
    f.close()

def remove_ion(fn, out):
    """
TER    1025        U A  47                                                      
HETATM 1026 MG    MG A 101      42.664  34.395  50.249  1.00 70.99          MG  
HETATM 1027 MG    MG A 201      47.865  33.919  48.090  1.00 67.09          MG 
    """
    pdb_fn = fn
    f = open(pdb_fn)

    ntxt = ''

    for l in f:
        element = l[76:78].strip().upper()
        if element in ['NA', 'MG']:
            continue
        else:
            ntxt += l
    f.close()

    f = open(out, 'w')
    f.write(ntxt)
    f.close()

def remove_water(fn, out):
    pdb_fn = fn
    f = open(pdb_fn)

    ntxt = ''

    for l in f:
        if l[17:20] in ['WAT', 'HOH']:
            continue
        else:
            ntxt += l
    f.close()

    f = open(out, 'w')
    f.write(ntxt)
    f.close()


def fix_rresnumes(fn,out):
    c = 1
    ntxt = ''
    i = open(fn)
    txt = i.read().split('\n')
    for l in txt:
        if l.startswith('ATOM') or l.startswith('HETATM') :
            #print l
            nl = l.replace( 'RA5', '  A') # RA should be the last!!!!
            nl = nl.replace('RA3', '  A')
            nl = nl.replace(' RA', '  A')
            nl = nl.replace(' rA', '  A')

            nl = nl.replace( 'RC5', '  C')
            nl = nl.replace('RC3', '  C')
            nl = nl.replace(' RC', '  C')
            nl = nl.replace(' rC', '  C')

            nl = nl.replace( 'RG5', '  G')
            nl = nl.replace('RG3', '  G')
            nl = nl.replace(' RG', '  G')
            nl = nl.replace(' rG', '  G')

            nl = nl.replace( 'RU5', '  U')
            nl = nl.replace('RU3', '  U')
            nl = nl.replace(' RU', '  U')
            nl = nl.replace(' rU', '  U')

            nl = nl.replace( 'RT5', '  T')
            nl = nl.replace('RT3', '  T')
            nl = nl.replace(' RT', '  T')
            nl = nl.replace(' rT', '  T')

            ntxt += nl + '\n'

        if l.startswith("END") or l.startswith("TER"):
            ntxt += l + '\n'

    i.close()
    o = open(out, 'w')
    o.write(ntxt)
    o.close()

def start():pass

if '__main__' == __name__:
    fn = 'test_data/image'
    print 'fn:', fn
    na = get_atom_lines(fn)
    print 'has_atom_line:', has_atom_line(na)

    fn = 'test_data/na.pdb'
    na = get_atom_lines(fn)
    print has_atom_line(na)
    res = get_all_res(na)
    print 'what is?', what_is(res)
    
    c = check_no_lines(na)
    print c
    res = get_all_res(na)
    #print res
    print 'non standard:', check_res_if_std_na(res)
    print 'is protein:', what_is(res)

    fn = 'test_data/prot.pdb'
    prot = open(fn).read().split('\n')
    res = get_all_res(prot)
    print 'non standard:', check_res_if_std_prot(res)   
    print 'is protein:', what_is(res)


    fn = 'test_data/rna-ru.pdb'
    prot = open(fn).read().split('\n')
    res = get_all_res(prot)
    print 'non standard:', check_res_if_supid_rna(res)
    print 'is protein:', what_is(res)

    fn = 'test_data/na_highAtomNum.pdb'
    renum_atoms(fn, '/tmp/out.pdb')

    fn = 'test_data/na_solvet_old_format.pdb'
    fix_op_atoms(fn, 'out.pdb')

    fn = 'test_data/na_solvet_old_format.pdb'
    fix_rresnumes(fn, 'out2.pdb')

    fn = 'test_data/na_solvet_old_format__.pdb'
    fix_rresnumes(fn, 'out2.pdb')
    remove_hydrogen('out2.pdb', 'out3.pdb')
    remove_ion('out3.pdb', 'out4.pdb')
    remove_water('out4.pdb', 'out5.pdb')
    renum_atoms('out5.pdb', 'out6.pdb')
    fix_op_atoms('out6.pdb', 'out7_.pdb')


    fn = 'test_data/1xjr.pdb'
    fix_rresnumes(fn, 'out2.pdb')
    remove_hydrogen('out2.pdb', 'out3.pdb')
    remove_ion('out3.pdb', 'out4.pdb')
    remove_water('out4.pdb', 'out5.pdb')
    renum_atoms('out5.pdb', 'out6.pdb')
    fix_op_atoms('out6.pdb', 'out7__.pdb')

    fn = 'test_data/decoy0165_amb.pdb'
    fix_rresnumes(fn, 'out2.pdb')
    remove_hydrogen('out2.pdb', 'out3.pdb')
    remove_ion('out3.pdb', 'out4.pdb')
    remove_water('out4.pdb', 'out5.pdb')
    renum_atoms('out5.pdb', 'out6.pdb')
    fix_O_in_UC('out6.pdb', 'out6.pdb')
    fix_op_atoms('out6.pdb', 'tmpout/decoy0165_amb_clx.pdb')

    fn = 'test_data/farna.pdb'
    fix_rresnumes(fn, 'out2.pdb')
    remove_hydrogen('out2.pdb', 'out3.pdb')
    remove_ion('out3.pdb', 'out4.pdb')
    remove_water('out4.pdb', 'out5.pdb')
    fix_op_atoms('out5.pdb', 'out7.pdb')
    renum_atoms('out7.pdb', 'out8farna.pdb')
