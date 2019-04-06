#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""RNAkb (previous Gromacs) utils.

A module with different functions needed for Gromacs/RNAkb merriage.

Marcin Magnus
Albert Bogdanowicz

(1) prepare groups and then (2) mdp score file.
"""

import re
import os

LIB_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + os.sep
VERBOSE = False
MDP_TEMPLATE = 'data/score_rnakb_orig.mdp'

GROMACS_ALLOWED_5 = ("H5T", "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'",
        "O4'", "C1'", "H1'", "N9", "C8", "H8", "N7", "C5", "C6", "N6", "H61",
        "H62", "N1", "C2", "H2", "N3", "C4", "C3'", "H3'", "C2'", "H2'1",
        "O2'", "HO'2", "O3'",
        # new atoms
        "O2", "O4", "O6", "N2", "N4",
        'OP1', 'OP2', 'OP3', 'H21', 'H22', "H2'", "H5'", "HO5'", "H5''", "HO2'",
        'P', 'H3', 'H1', 'H6', 'H5', 'H42', 'H41',
        )

GROMACS_ALLOWED_MIDDLE = ("P", "O1P", "O2P", "O5'", "C5'", "H5'1", "H5'2",
       "C4'", "H4'", "O4'", "C1'", "H1'", "N9", "C8", "H8", "N7", "C5", "C6",
        "N6", "H61", "H62", "N1", "C2", "H2", "N3", "C4", "C3'", "H3'", "C2'",
        "H2'1", "O2'", "HO'2", "O3'",
        # new atoms
        "O2", "O4", "O6", "N2", "N4",
        'OP1', 'OP2', 'OP3', 'H21', 'H22', 'O6', "H2'", "H5'", "HO5'", "H5''", "HO2'",
        'P', 'H3', 'H1', 'H6', 'H5', 'H42', 'H41',
        )

GROMACS_REQUIRED_GROUPS = ("aP", "aC4s", "aC2", "aC4", "aC6",
                           "uP", "uC4s", "uC2", "uC4", "uC6",
                           "gP", "gC4s", "gC2", "gC4", "gC6",
                           "cP", "cC4s", "cC2", "cC4", "cC6",
                           "RNA_5pt", "other"
        )


def get_res_num(line):
    """Extract residue number from a line of PDB file

    Arguments:
      * line = ATOM line from a PDB file

    Output:
      * residue number as an integer
    """
    return int(''.join([x for x in line[22:27] if x.isdigit()]))


def get_res_code(line):
    """Get residue code from a line of a PDB file
    """
    return line[17:20]

def set_res_code(line, code):
    """Set residue code from a line of a PDB file
    """
    return line[:17] + code.rjust(3) + line[20:]

def make_rna_gromacs_ready(pdb_string, verbose=VERBOSE):
    """GROMACS has some special requirements for PDB files.

    Arguments:
      * pdb_string = contents of PDB file as a string

    Output:
      * new PDB returned as a string

    (!!!) # hmm... [ RA5 ] will not be detected based on it (!?)
    Hmm.. becase it detects if the structure is already prepared.
    """
    #pdb_string = resname_3to1(pdb_string)
    #pdb_string = remove_hetatms(pdb_string)
    result = []
    pdb_lines = pdb_string.split('\n')

    # find smallest residue number
    min_res = min(list(map(get_res_num,
        [l for l in pdb_lines if l.startswith('ATOM')])))
    max_res = max(list(map(get_res_num,
        [l for l in pdb_lines if l.startswith('ATOM')])))


    for l in pdb_lines:
        if l.startswith('ATOM') and l[19] in ('A', 'U', 'C', 'G'): # hmm... [ RA5 ] will not be detected based on it (!?)
            res = get_res_code(l)        
            if res.startswith('R') and res.endswith('5') : # it's RX5 file so skip fixing
                if verbose:
                    print('-- already gromacs ready')
                return pdb_string

            l = l.replace('*', '\'')
            l = l.replace('O1P', 'OP1')
            l = l.replace('O2P', 'OP2')

            res_num = get_res_num(l)
            atom_type = l.split()[2].strip()
            # remove P OP1 OP2
            if res_num == min_res and atom_type == 'P':
                continue
            if res_num == min_res and atom_type == 'OP1':
                continue
            if res_num == min_res and atom_type == 'OP2':
                continue

            # convert G -> RG5, RG3
            if res_num == min_res: # RG5
                l = set_res_code(l, 'R' + get_res_code(l).strip() + '5')
            elif res_num == max_res: # RG3
                l = set_res_code(l, 'R' + get_res_code(l).strip() + '3')
            else:
                l = set_res_code(l, ' R' + get_res_code(l).strip())

            if res_num == min_res:
                if atom_type in GROMACS_ALLOWED_5:
                    result.append(l)
                else:
                    print(('Wrong start line: ', l, atom_type))
            else:
                if atom_type in GROMACS_ALLOWED_MIDDLE:
                    result.append(l)
                else:
                    print(('Wrong middle line: ', l, atom_type))
        else:
            result.append(l)
    return '\n'.join(result)

def make_rna_rnakb_ready(pdb_string, verbose=VERBOSE):
    """RNAkb read (difference between this function and 
    make_rna_gromacs_ready is ignoring R5U etc. RNAkb does not treat
    them differently so there is no point to distinguish them.

    Arguments:
      * pdb_string = contents of PDB file as a string

    Output:
      * new PDB returned as a string
    """
    #pdb_string = resname_3to1(pdb_string)
    #pdb_string = remove_hetatms(pdb_string)
    result = []
    pdb_lines = pdb_string.split('\n')

    # find smallest residue number
    min_res = min(list(map(get_res_num,
        [l for l in pdb_lines if l.startswith('ATOM')])))
    max_res = max(list(map(get_res_num,
        [l for l in pdb_lines if l.startswith('ATOM')])))

    for l in pdb_lines:
        if l.startswith('ATOM'):# and l[19] in ('A', 'U', 'C', 'G'):
            res = get_res_code(l)
            #if res.startswith('R') and res.endswith('5') : # it's RX5 file so skip fixing
            #    if verbose:
            #        print '-- already gromacs ready'
            #    return pdb_string

            l = l.replace('*', '\'')
            l = l.replace('O1P', 'OP1')
            l = l.replace('O2P', 'OP2')

            res_num = get_res_num(l)
            atom_type = l.split()[2].strip()
            # remove P OP1 OP2
            if res_num == min_res and atom_type == 'P':
                continue
            if res_num == min_res and atom_type == 'OP1':
                continue
            if res_num == min_res and atom_type == 'OP2':
                continue

            # convert G -> RG5, RG3
            #if res_num == min_res: # RG5
            #    l = set_res_code(l, 'R' + get_res_code(l).strip() + '5')
            #elif res_num == max_res: # RG3
            #    l = set_res_code(l, 'R' + get_res_code(l).strip() + '3')
            #else:
            l = set_res_code(l, ' R' + get_res_code(l).strip().replace('R','').replace('3','').replace('5','')) # 
            if res_num == min_res:
                if atom_type in GROMACS_ALLOWED_5:
                    result.append(l)
                else:
                    print(('Wrong start line: ', l, atom_type))
            else:
                if atom_type in GROMACS_ALLOWED_MIDDLE:
                    result.append(l)
                else:
                    print(('Wrong middle line: ', l, atom_type))
        else: # keep TER, etc.
            result.append(l)
    return '\n'.join(result)

def fix_gromacs_gro(path, verbose=False):
    """It's probably a bug in GROMACS, but box coordinates in gro files are
    not always separated by spaces. This function guesses how it should be
    separated and inserts spaces.

    Arguments:
      * path = path to gro file

    Output:
      * file is overwritten with a corrected one
    """
    f = open(path)
    gro = f.read()
    f.close()
    gro_lines = gro.split('\n')
    last_line = gro_lines[-2]

    # check if there are a space
    if last_line.find(' ') == -1:
        dots = [i.start() for i in re.finditer('\\.', last_line)]
        # next 4 lines are a guess, I hope it works
        digits = len(last_line[dots[2]:])
        box = [last_line[:dots[0] + digits],
            last_line[dots[0] + digits:dots[1] + digits], last_line[dots[1] + digits:]]
        gro_lines = gro_lines[:-2]
        gro_lines.append(' '.join(box))
        gro_lines.append('')
        f = open(path, 'w')
        f.write('\n'.join(gro_lines))
        f.close()


def fix_gromacs_ndx(path):
    """Sometimes, GROMACS index has some atoms in more than one group, or
    doesn't have all the groups grompp requires. This function fixes that.

    Arguments:
      * path = path to index file

    Output:
      * index is overwritten with a corrected one
    """
    f = open(path)
    index = f.read()
    f.close()
    # split file into system group and the rest
    system_group = index.split('[ System ]')[1].split('[')[0]
    other_groups = ['[' + i for i in
            index.split('[ System ]')[1].split('[')[1:]]
    # remove duplicate numbers
    taken_atoms = []
    for g in range(len(other_groups)):
        header = other_groups[g].split('\n')[0]
        group = other_groups[g].split('\n')[1].split()
        group = [a for a in group if a not in taken_atoms]
        taken_atoms.extend(group)
        other_groups[g] = header + '\n' + ' '.join(group) + '\n'
    # build result, part 1
    result = ['[ System ]' + system_group]
    result.extend(other_groups)
    # add missing groups, leave them empty
    headers = [g.split('\n')[0][2:-2] for g in other_groups]
    missing_headers = [h for h in GROMACS_REQUIRED_GROUPS if h not in headers]
    result.extend(['[ %s ]\n' % h for h in missing_headers])
    # write result to file
    result = ''.join(result)
    f = open(path, 'w')
    f.write(result)
    f.close()


def prepare_groups(fn, gr_fn, potential='aa', verbose=False):
    """ Prepare an index for fn file. gr_fn is a file where gtxt is saved in.

    Get seq and uniq & sort it.
    ``['RG5', 'RA', 'RA', 'RA', 'RG', 'RU', 'RA', 'RA', 'RC3'] set(['RU', 'RG', 'RC3', 'RG5', 'RA'])``

    @todo RG5 etc -- done!

    gtxt::

     del 1
     r RU & a C1'
     name 1 uC1s
     r RU & a C2
     name 2 uC2
     r RU & a C2'
     name 3 uC2s
     ...

    return, gtxt (groups_txt), energygrps . The result is saved to g_fn.
    energygrps:
    ['uP', 'uC4s', 'uC2', 'uC4', 'uC6', 'gP', 'gC4s', 'gC2', 'gC4', 'gC6', 'aP', 'aC4s', 'aC2', 'aC4', 'aC6']
    gtxt:
    RA
    del 1
    r RU & a P
    name 1 uP
    r RU & a C4s
    name 2 uC4s
    r RU & a C2
    name 3 uC2
    r RU & a C4
    [..]
    r RA & a C6
    name 15 aC6
    1|2|3|4|5|6|7|8|9|10|11|12|13|14|15
    name 16 RNA_5pt
    0 & ! 16
    name 17 other
    q
    """
    p = pf.PDBFile(pdb_path = fn) 
    seq = p.seq_from_amber_like_pdb().split()
    seq_uniq_sorted = set(seq)
    if verbose: print(('seq:', seq, seq_uniq_sorted))
    seq_rnakb_order = []
    if 'RA' in seq_uniq_sorted:
        seq_rnakb_order.append('RA')
    if 'RU' in seq_uniq_sorted:
        seq_rnakb_order.append('RU')
    if 'RG' in seq_uniq_sorted:
        seq_rnakb_order.append('RG')
    if 'RC' in seq_uniq_sorted:
        seq_rnakb_order.append('RC')
    if verbose:print(('seq_rnakb_order', seq_rnakb_order))
    
    gtxt = 'del 1\n'
    c = 1

    if potential == 'aa':
        # rg
        rg_atoms = "C3',C5,C4,C6,C8,O2',P,C2',O5',C5',C1',O3',O6,N2,N3,N1,N7,N9,C2,C4',O4',O2P,O1P".split(',')
        rg_atoms2 = ['g' + a.strip().replace("'", 's') for a in rg_atoms]
        # ru
        ug_atoms =  "C1',C2,C2',C3',C4,C4',C5,C5',C6,N1,N3,O2,O2',O3',O4,O4',O5',O1P,O2P,P".split(",")
        ug_atoms2 = ['u' + a.strip().replace("'", 's') for a in ug_atoms]
        # ag
        ag_atoms = "O3',O2',N7,N1,N3,N9,C2',O5',N6,C5',C1',C2,C6,C5,C4,O4',C4',C8,C3',P,O1P,O2P".split(',')
        ag_atoms2 = ['a' + a.strip().replace("'", 's') for a in ag_atoms]
        # cg
        cg_atoms = "C2',O2',O2P,O1P,C5',O5',C4,O2,C3',C2,O3',N4,N3,N1,P,C1',O4',C4',C5,C6".split(',')
        cg_atoms2 = ['c' + a.strip().replace("'", 's') for a in cg_atoms]
        #N1, N3, O4', C5', O3', C2', C4, C1', O5', O1P, C4', C6, C5, C2, C3', P, O2P, O2, O4, O2'
        if verbose: print(('len-s:', len(rg_atoms), len(cg_atoms), len(ag_atoms), len(ug_atoms)))
    elif potential == '5pt':
        # rg
        rg_atoms = ["P", "C4'", "C2", "C4", "C6"]
        ug_atoms = ["P", "C4'", "C2", "C4", "C6"]
        ag_atoms = ["P", "C4'", "C2", "C4", "C6"]                
        cg_atoms = ["P", "C4'", "C2", "C4", "C6"]                
        rg_atoms2 = ['g' + a.strip().replace("'", 's') for a in rg_atoms]
        ug_atoms2 = ['u' + a.strip().replace("'", 's') for a in ug_atoms]
        ag_atoms2 = ['a' + a.strip().replace("'", 's') for a in ag_atoms]
        cg_atoms2 = ['c' + a.strip().replace("'", 's') for a in cg_atoms]
    else:
        raise Exception('Unknown potential, use all or 5pt')

    energygrps = []

    for r in seq_rnakb_order:
        if r == 'RA':
            for x, y in zip(ag_atoms, ag_atoms2):
                gtxt += 'r RA & a %s\n' % x
                gtxt += 'name %i %s\n' % (c, y)
                c += 1
            energygrps.extend(ag_atoms2)

        if r == 'RU':
            for x, y in zip(ug_atoms, ug_atoms2):
                gtxt += 'r RU & a %s\n' % x
                gtxt += 'name %i %s\n' % (c, y)
                c += 1
            energygrps.extend(ug_atoms2)

        if r == 'RG':
            for x, y in zip(rg_atoms, rg_atoms2):
                gtxt += 'r RG & a %s\n' % x
                gtxt += 'name %i %s\n' % (c, y)
                c += 1
            energygrps.extend(rg_atoms2)

        if r == 'RC':
            for x, y in zip(cg_atoms, cg_atoms2):
                gtxt += 'r RC & a %s\n' % x
                gtxt += 'name %i %s\n' % (c, y)
                c += 1
            energygrps.extend(cg_atoms2)

    gtxt += '|'.join([str(x) for x in range(1,c)])
    gtxt += '\nname %i RNA_%s' % (c, potential) # @todo
    gtxt += '\n0 & ! %i' % (c)
    gtxt += '\nname %i other' % (c + 1)
    gtxt += '\nq\n'
    if verbose: print(gtxt)

    with open(gr_fn, 'w') as f:
        f.write(gtxt)

    if verbose: print(('energygrps', energygrps))
    return gtxt, energygrps, seq_rnakb_order


def format_score_mdp(mdp_out, energygrps, seq, verbose=False):
    """Get a template score mdp and replace energygrps
    (it can be generated with prepare_groups)
    and energygrp_table
    """
    # load template
    with open(LIB_PATH + 'rnakb_utils/' + MDP_TEMPLATE, 'r') as f:
        txt = f.readlines()

    with open(LIB_PATH + 'rnakb_utils/data/rnakb_all.txt', 'r') as f:
        pairs = [i.strip() for i in f.readlines()]

    nmdp = ''
    for l in txt:
        if l.startswith('energygrps'):
            l = 'energygrps               = ' + ' '.join(energygrps) + ' other'
            nmdp += l
        elif l.startswith('energygrp_table'):
            d = ''
            for x in energygrps:  # ugly :-(
                for y in energygrps:
                    s = '%s_%s' % (x, y) # Library file table_uP_aP.xvg not found in current dir nor in your GMXLIB path.
                    s_reverse = '%s_%s' % (y, x) # try: /home/mqapRNA/mqaprna_env/db/RNA_5pt_full_sc1/table_aP_uP.xvg
                    if s in pairs:
                        d += '%s ' % s.replace('_', ' ')  # '_' -> ' '
            l = 'energygrp_table          = ' + d.strip()
            nmdp += l
        elif l.startswith('energygrp_excl'):
            l = 'energygrp_excl           =  other other ' + ' other '.join(energygrps) + ' other'
            nmdp += l
        else:
            nmdp += l
    if verbose: print(nmdp)

    with open(mdp_out, 'w') as f:
        f.write(nmdp)

# main
if __name__ == '__main__':
    #fn = 'test_data/decoy0165_amb_clx.pdb'
    #fn = 'test_data/1duq.pdb'
    fn = 'test_data/cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_noC.pdb'
    pdblines = make_rna_gromacs_ready(open(fn).read())
    print(pdblines)
    with open('test_output/gromacs_ready.pdb', 'w') as f:
        f.write(pdblines)

    fn = 'test_data/cat_chunk003_2r8s_5kcycles_decoys_nonativefrags.cluster1.0_clean_noC.pdb'
    pdblines = make_rna_rnakb_ready(open(fn).read())
    fready = 'test_output/rnakb_ready.pdb'
    print(pdblines)
    with open(fready, 'w') as f:
        f.write(pdblines)

    # prepare groups
    groups_txt, energygrps, seq_uniq = prepare_groups(fready, LIB_PATH + '/rnakb_utils/test_output/groups.txt', potential='5pt', verbose=True)
    print(groups_txt)
    fout = 'test_data/out.mdp'
