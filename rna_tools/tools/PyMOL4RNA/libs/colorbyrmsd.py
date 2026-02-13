'''
http://pymolwiki.org/index.php/ColorByRMSD

Original Authors: Shivender Shandilya; Jason Vertrees
Complete rewrite by Thomas Holder
Adaption for RNA by Marcin Magnus

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def colorbyrmsd(mobile, target, doAlign=1, doPretty=1, guide=1, method='align', quiet=1):
    '''
    ! super does not work with RNA, so we use align instead, which is more robust for RNA !
    
    DESCRIPTION

    Align two structures and show the structural deviations in color to more
    easily see variable regions.

    Colors each mobile/target atom-pair by distance (the name is a bit
    misleading).

    Modifies the B-factor columns in your original structures.

    This script allows you to color two structures by Root Mean Square Deviation (RMSD). The distances between aligned C-alpha atom pairs are stored as B-factors of these residues, which are colored by a color spectrum, with blue specifying the minimum pairwise RMSD and red indicating the maximum. Unaligned residues are colored gray.

    To change ranges of spectrum simply run the spectrum command again with minimum and maximum arguments.

    PyMOL> spectrum b, blue_red, b > -0.5, minimum=0, maximum=3.0

    https://pymol-users.narkive.com/auBmoGg9/pymol-colorbyrmsd

ARGUMENTS

    mobile = string: atom selection for mobile atoms
    
    target = string: atom selection for target atoms

    doAlign = 0 or 1: Superpose selections before calculating distances
    {default: 1}

    doPretty = 0 or 1: Show nice representation and colors {default: 1}

EXAMPLE

    fetch 1ake 4ake, async=0
    remove chain B
    colorbyrmsd 1ake, 4ake
    '''
    from chempy import cpv
    ic = None
    try:
        from icecream import ic as _ic
        import sys
        _ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
        _ic.configureOutput(prefix='')
        ic = _ic
    except ImportError:
        pass

    doAlign, doPretty = int(doAlign), int(doPretty)
    guide, quiet = int(guide), int(quiet)
    aln, seleboth = '_aln', '_objSelBoth'

    mobile_full = '(%s)' % mobile
    target_full = '(%s)' % target

    try:
        align = cmd.keyword[method][0]
    except:
        print(' Error: no such method:', method)
        raise CmdException

    def apply_guide_atoms(selection):
        if not guide:
            return selection
        # try PyMOL's default guide atoms (CA for proteins, etc.) first
        candidates = [
            ('guide', '(%s) and guide' % selection),
            # fall back to nucleic-acid specific atoms
            ("RNA C4*", '(%s) and polymer.nucleic and name C4*' % selection),
            ('RNA P', '(%s) and polymer.nucleic and name P' % selection),
        ]
        for label, candidate in candidates:
            if cmd.count_atoms(candidate):
                if not quiet:
                    print(' ColorByRMSD: using %s atoms for %s' % (label, selection))
                return candidate
        if not quiet:
            print(' ColorByRMSD: no guide atoms found for %s, using full selection' % selection)
        return '(%s)' % selection

    mobile_sel, target_sel = mobile_full, target_full
    if guide:
        mobile_sel = apply_guide_atoms(mobile_full)
        target_sel = apply_guide_atoms(target_full)

    try:
        if doAlign:
            # superpose
            align(mobile_sel, target_sel)

        # get alignment without superposing
        align(mobile_sel, target_sel, cycles=0, transform=0, object=aln)
    except:
        print(' Error: Alignment with method %s failed' % (method))
        raise CmdException

    raw_aln = cmd.get_raw_alignment(aln)
    for idx1, idx2 in raw_aln:
        print('%s`%d -> %s`%d' % tuple(idx1 + idx2))

    if cmd.select(seleboth, '(%s) or (%s)' % (mobile_full, target_full)) == 0:
        raise CmdException('Selections have no atoms: %s or %s' % (mobile_full, target_full))

    # gather coordinates for atoms referenced in the alignment
    idx2coords = dict()
    idx2residue = dict()
    aligned_atoms = set()
    for col in raw_aln:
        aligned_atoms.update(col)

    missing_atoms = []
    for model, index in aligned_atoms:
        sele = '(%s and index %d)' % (model, index)
        if cmd.count_atoms(sele) == 0:
            missing_atoms.append('%s`%d' % (model, index))
            continue
        idx2coords[(model, index)] = cmd.get_atom_coords(sele, state=-1)
        atom_model = cmd.get_model(sele)
        if atom_model.atom:
            atom = atom_model.atom[0]
            idx2residue[(model, index)] = (atom.segi, atom.chain, atom.resi)

    if missing_atoms and not quiet:
        print(' ColorByRMSD: skipping atoms with missing coordinates:', ', '.join(missing_atoms))

    if cmd.count_atoms('?' + aln, 1, 1) == 0:
        # this should ensure that "aln" will be available as selectable object
        cmd.refresh()

    b_dict = dict()
    for col in raw_aln:
        if ic:
            ic(col)
        assert len(col) == 2
        if col[0] not in idx2coords or col[1] not in idx2coords:
            if not quiet:
                print(' ColorByRMSD: alignment pair skipped due to missing coords: %s`%d - %s`%d' % tuple(col[0] + col[1]))
            continue
        b = cpv.distance(idx2coords[col[0]], idx2coords[col[1]])
        if ic:
            ic(b)
        for idx in col:
            b_dict[idx] = b

    if not b_dict:
        raise CmdException('No aligned atom pairs with coordinates found')

    residue_b = dict()
    for idx, b in b_dict.items():
        residue = idx2residue.get(idx)
        if not residue:
            continue
        key = (idx[0],) + residue  # (model, segi, chain, resi)
        residue_b[key] = b

    if not residue_b:
        raise CmdException('No residue mapping available for aligned atoms')

    cmd.alter(
        seleboth,
        'b = residue_b.get((model, segi, chain, resi), -1)',
        space={'residue_b': residue_b},
    )

    if cmd.count_atoms(seleboth) == 0:
        raise CmdException('Temporary selection %s is empty' % seleboth)

    if doPretty:
        cmd.orient(seleboth)
        cmd.show_as('cartoon', 'byobj ' + seleboth)
        cmd.color('gray', seleboth)
        cmd.spectrum('b', 'blue_red', seleboth + ' and b > -0.5')

    if ic:
        ic(b_dict.values())
    if not quiet:
        print(" ColorByRMSD: Minimum Distance: %.2f" % (min(b_dict.values())))
        print(" ColorByRMSD: Maximum Distance: %.2f" % (max(b_dict.values())))
        print(" ColorByRMSD: Average Distance: %.2f" % (sum(b_dict.values()) / len(b_dict)))

    cmd.delete(aln)
    cmd.delete(seleboth)

cmd.extend('colorbyrmsd', colorbyrmsd)

# tab-completion of arguments
cmd.auto_arg[0]['colorbyrmsd'] = cmd.auto_arg[0]['align']
cmd.auto_arg[1]['colorbyrmsd'] = cmd.auto_arg[1]['align']

# vi: ts=4:sw=4:smarttab:expandtab
