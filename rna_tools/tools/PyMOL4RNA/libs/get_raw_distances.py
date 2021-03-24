"""
Authors: (c) 2012 Takanori Nakane and Thomas Holder License: BSD-2-Clause

Modified: Marcin Magnus 2020

Source <http://pymolwiki.org/index.php/get_raw_distances>
"""

from pymol import cmd, CmdException

def aa3to1(aaa):
    """based on https://pymolwiki.org/index.php/Aa_codes"""
    if len(aaa) != 3:  # aaa is 'G', like for RNA ;-)
        return '' # dont do it for rna test aaa
    one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
                 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
                 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
                 'GLY':'G', 'PRO':'P', 'CYS':'C'}
    return one_letter[aaa]


def aa1to3(a):
    """based on https://pymolwiki.org/index.php/Aa_codes"""
    three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
    'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
    'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
    'G':'GLY', 'P':'PRO', 'C':'CYS'}
    return three_letter[a]


def unid(object, index, format="csv"):
    """
    l = [['MET', '1'], ['MET', '1'], ['MET', '1'], ['MET', '1'], ['MET', '1'], ['THR', '2'], ['MET', '1'], ['THR', '2'], ['MET', '1'], ['THR', '2'], ['MET', '1'], ['THR', '2'], ['THR', '2'], ['THR', '2'], ['THR', '2'], ['THR', '2'], ['THR', '2'], ['THR', '2'], ['THR', '2'], ['THR', '2'], ['MET', '1'], ['THR', '2'], ['MET', '1'], ['THR', '2'], ['THR', '2'], ['THR', '3'], ['THR', '2'], ['THR', '3'], ['THR', '2'], ['THR', '3'], ['THR', '2'], ['THR', '3'], ['THR', '3'], ['THR', '3'], ['THR', '3'], ['THR', '3'], ['THR', '2'], ['THR', '3'], ['THR', '2'], ['THR', '3'], ['THR', '3'], ['THR', '3'], ['THR', '3'], ['THR', '3'], ['THR', '3'], ['SER', '4'], ['THR', '3'], ['SER', '4'], ['THR', '3'], ['SER', '4'], ['THR', '3'], ['SER', '4'], ['THR', '3'], ['SER', '4'], ['THR', '3'], ['SER', '4'], ['THR', '2'], ['G', '52'], ['THR', '2'], ['G', '52'], ['THR', '2'], ['G', '52'], ['THR', '2'], ['G', '52'], ['THR', '2'], ['A', '53'], ['THR', '2'], ['A', '53'], ['THR', '2'], ['A', '53'], ['THR', '2'], ['A', '53'], ['THR', '2'], ['A', '53'], ['THR', '2'], ['A', '53'], ['THR', '3'], ['A', '53'], ['THR', '3'], ['A', '53'], ['THR', '3'], ['A', '53'], ['THR', '3'], ['A', '53'], ['THR', '3'], ['A', '53'], ['THR', '3'], ['A', '53'], ['THR', '3'], ['U', '54']]

"""
    # ('yC_5lj3_U6', 1413)
    #for r in
    # print(cmd.iterate(, 'print(resi)'))
    #cmd.do("python('l = [];')")
    selection = object + ' and index ' + str(index)
    from pymol import stored
    stored.l = []
    #cmd.iterate(selection, 'l.append([resn, resi])')
    cmd.iterate(selection, 'stored.l.append([resn, resi])')
    #print(stored.l)
    #print(l)
    #l = cmd.do('l')
    
    if format == 'csv': 
        mol = object.split('_')[0] # yC_5lj3_Prp8 ## TODO
        if stored.l:
            # ['THR', '3'], ['A', '53'], 
            resn = stored.l[-1][0]
            resi = stored.l[-1][1]
            if len(resn) == 3: # process aa from THR to T
                resn = aa3to1(resn) # else keep it C, G, U, A, etc
            x = mol + '-' + resn + str(resi)
        else:
            x = ''
    else:
        x = object + ' and resi ' +  str(aa3to1(stored.l[-1][0])) + str(stored.l[-1][1])
        # print(x)
    return (x)
    #: # 'id ' + str(1411):
    #    print(r)


def get_raw_distances(names='', state=1, selection='all', quiet=1, filename='intrs.csv'):
    '''
DESCRIPTION

    Get the list of pair items from distance objects. Each list item is a
    tuple of (index1, index2, distance).

    Based on a script from Takanori Nakane, posted on pymol-users mailing list.
    http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10143.html

ARGUMENTS

    names = string: names of distance objects (no wildcards!) {default: all
    measurement objects}

    state = integer: object state {default: 1}

    selection = string: atom selection {default: all}

    quiet = boolen

    filename = if this is '' then dont save any file at all, be default 'intrs.csv'

SEE ALSO

    select_distances, cmd.find_pairs, cmd.get_raw_alignment
    '''
    #foo = cmd.do('l = [];') ## ugly hack!
    from chempy import cpv

    state, quiet = int(state), int(quiet)
    if state < 1:
        state = cmd.get_state()

    valid_names = cmd.get_names_of_type('object:measurement')
    if names == '':
        names = ' '.join(valid_names)
    else:
        for name in names.split():
            if name not in valid_names:
                print(' Error: no such distance object: ' + name)
                raise CmdException

    raw_objects = cmd.get_session(names, 1, 1, 0, 0)['names']

    xyz2idx = {}
    cmd.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model,index)',
                      space=locals())

    r = []
    rresi = []
    for obj in raw_objects:
        try:
            points = obj[5][2][state - 1][1]
            if not quiet:
                print(points)
            if points is None:
                raise ValueError
        except (KeyError, ValueError):
            continue
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i + 3])
            xyz2 = tuple(points[i + 3:i + 6])
            try:
                r.append((xyz2idx[xyz1], xyz2idx[xyz2], cpv.distance(xyz1, xyz2)))
                # (('yC_5lj3_U6', 1183)
                if not quiet:
                    print(' get_raw_distances: ' + str(r[-1]))
                rresi.append([unid(xyz2idx[xyz1][0], xyz2idx[xyz1][1]), unid(xyz2idx[xyz2][0], xyz2idx[xyz2][1])])
                if not quiet:
                    print('  ', unid(xyz2idx[xyz1][0], xyz2idx[xyz1][1]), '<->', unid(xyz2idx[xyz2][0], xyz2idx[xyz2][1]))
            except KeyError:
                if quiet < 0:
                    print(' Debug: no index for %s %s' % (xyz1, xyz2))

    if rresi:
        if filename: # if empty filename then don't save it
            with open(filename, 'w') as f:
                for r in rresi:
                    f.write(r[0] + ',' + r[1] + '\n')
            print('File saved:', filename)
    return r, rresi


def select_distances(names='', name='sele', state=1, selection='all', cutoff=-1, quiet=1):
    '''
DESCRIPTION

    Turns a distance object into a named atom selection.

ARGUMENTS

    names = string: names of distance objects (no wildcards!) {default: all
    measurement objects}

    name = a unique name for the selection {default: sele}

SEE ALSO

    get_raw_distances
    '''
    state, cutoff, quiet = int(state), float(cutoff), int(quiet)

    sele_dict = {}
    distances = get_raw_distances(names, state, selection)
    for idx1, idx2, dist in distances:
        if cutoff <= 0.0 or dist <= cutoff:
            sele_dict.setdefault(idx1[0], set()).add(idx1[1])
            sele_dict.setdefault(idx2[0], set()).add(idx2[1])

    cmd.select(name, 'none')
    tmp_name = cmd.get_unused_name('_')

    r = 0
    for model in sele_dict:
        cmd.select_list(tmp_name, model, list(sele_dict[model]), mode='index')
        r = cmd.select(name, tmp_name, merge=1)
        cmd.delete(tmp_name)

    if not quiet:
        print(' Selector: selection "%s" defined with %d atoms.' % (name, r))
    return r

cmd.extend('get_raw_distances', get_raw_distances)
cmd.extend('select_distances', select_distances)

_auto_arg0_distances = [
    lambda: cmd.Shortcut(cmd.get_names_of_type('object:measurement')),
    'distance object', '']

cmd.auto_arg[0].update([
    ('get_raw_distances', _auto_arg0_distances),
    ('select_distances', _auto_arg0_distances),
])
