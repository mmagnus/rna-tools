#!/usr/bin/env python
"""
Read for more interesting functions https://daslab.stanford.edu/site_data/docs_pymol_rhiju.pdf
"""
import tempfile
import math
import subprocess
import os
from itertools import izip

from rna_pdb_tools.rna_pdb_tools_lib import RNAStructure

try:
    RNA_PDB_TOOLS
    EXECUTABLE
except NameError:
    RNA_PDB_TOOLS = os.environ.get('RNA_PDB_TOOLS')
    EXECUTABLE="/bin/zsh"
    SOURCE=""

def exe(cmd, verbose=False):
    """Helper function to run cmd. Using in this Python module."""
    if verbose: print('cmd:' + cmd)
    o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         executable=EXECUTABLE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


def color_by_text(txt):
    """Helper function used for color-coding based on residue indexes ranges."""
    for t in txt.strip().split('\n'):
        color, resi = t.replace('color ', '').split(',')
        print((color, resi))
        cmd.color(color.strip(), resi.strip())


def rp():
    """Represent your RNA."""
    cmd.hide("sticks", "all")
    cmd.hide("lines", "all")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 2)
    cmd.set("cartoon_ladder_mode", 1)


def rs():
    """    The function creates super-cool cartoon-like RNA and colors each structure as a rainbow.
    Good to view aligned structures in a grid.

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/rs.png
    """
    cmd.hide("sticks", "all")
    cmd.hide("lines", "all")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 2)
    cmd.set("cartoon_ladder_mode", 2)
    cmd.set("cartoon_ring_transparency", 0.30)
    cmd.spectrum()

    obj_list = cmd.get_names('objects')

    colours = ['rainbow']
    ncolours = len(colours)
    # Loop over objects
    i = 0
    for obj in obj_list:
        print "  ", obj, colours[i]
        cmd.spectrum('count', colours[i], obj)
        i = i+1
        if(i == ncolours):
           i = 0


def rcomp():
    """RNA like in papers ;-)

    Similar to rc() but this time it colors each (and every) structure in different colour.
    Great on viewing-comparing superimposed structures.

    """
    cmd.hide("sticks", "all")
    cmd.hide("lines", "all")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 2)
    cmd.set("cartoon_ladder_mode", 2)
    cmd.set("cartoon_ring_transparency", 0.30)

    obj_list = cmd.get_names('objects')

    colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
           'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', \
           'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',    \
           'wheat', 'white', 'grey' ]
    ncolours = len(colours)

           # Loop over objects
    i = 0
    for obj in obj_list:
        print "  ", obj, colours[i]
        cmd.color(colours[i], obj)
        i = i+1
        if(i == ncolours):
           i = 0


def align_all( subset = [] ):
  """
  Superimpose all open models onto the first one.
  This may not work well with selections.

  This function is probably taken from https://daslab.stanford.edu/site_data/docs_pymol_rhiju.pdf
  """
  print """This returns a list with 7 items:

    RMSD after refinement
    Number of aligned atoms after refinement
    Number of refinement cycles
    RMSD before refinement
    Number of aligned atoms before refinement
    Raw alignment score
    Number of residues aligned """

  AllObj=cmd.get_names("all")
  for x in AllObj[1:]:
    #print(AllObj[0],x)
    subset_tag = ''
    if isinstance( subset, int ):
      subset_tag = ' and resi %d' % subset
    elif isinstance( subset, list ) and len( subset ) > 0:
      subset_tag = ' and resi %d' % (subset[0])
      for m in range( 1,len(subset)): subset_tag += '+%d' % subset[m]
    elif isinstance( subset, str ) and len( subset ) > 0:
      subset_tag = ' and %s' % subset
    values = cmd.align(x+subset_tag,AllObj[0]+subset_tag)
    print AllObj[0], x, ' '.join([str(v) for v in values]), '-- RMSD', values[3], ' of ', values[6], 'residues'
    cmd.zoom()


def get_pdb():
    """Get PDB content of selection.

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/pdb.png"""
    tmpfn = '/tmp/pymol_get_pdb.pdb'
    cmd.save(tmpfn, '(sele)')
    s = RNAStructure(tmpfn)
    for l in s.lines:
        print(l)


def clarna():
    """Get contacts classification of the selected fragment based on ClaRNA.

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/clarna.png
    """
    f = tempfile.NamedTemporaryFile(delete=False) # True)
    cmd.save(f.name + '.pdb', '(sele)')
    out, err = exe(SOURCE + " && " + CLARNA_RUN + " -ipdb " + f.name + '.pdb -bp+stack')
    print('\n'.join(out.split('\n')[1:]))  # to remove first line of py3dna /tmp/xxx
    if err:
        print(err)
    f.close()


def get_seq():
    """Get contacts classification based on ClaRNA.

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/ss.png
    """
    f = tempfile.NamedTemporaryFile(delete=False) # True)
    cmd.save(f.name, '(sele)')
    out, err = exe('source ~/.zshrc && ' + RNA_PDB_TOOLS + '/bin/rna_pdb_toolsx.py --get_seq ' + f.name)
    print(out)
    if err:
        print(err)
    f.close()

def ss():
    """Get Secondary Structure of (sele) based on py3dna.py.

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/ss.png
    """
    f = tempfile.NamedTemporaryFile(delete=False) # True)
    cmd.save(f.name, '(sele)')
    out, err = exe(RNA_PDB_TOOLS + '/bin/rna_x3dna.py ' + f.name)
    print('\n'.join(out.split('\n')[2:]))  # to remove first line of py3dna /tmp/xxx
    if err:
        print(err)
    f.close()


def ss_all():
    """The same as ss() but for all objects."""
    subset = "*"
    AllObj = cmd.get_names("all")
    # print AllObj
    for name in AllObj[:]:
        if not name.startswith('_align'):
            print('> ' + name)
            f = tempfile.NamedTemporaryFile(delete=False) # True)
            cmd.save(f.name, name)
            out, err = exe(RNA_PDB_TOOLS + '/bin/rna_x3dna.py ' + f.name)
            print('\n'.join(out.split('\n')[2:]))  # to remove first line of py3dna /tmp/xxx
            # hide this line: is >tmpGCszi7 nts=4 [tmpGCszi7] -- secondary structure derived by DSSR
            if err:
                print(err)
            f.close()
    print('-- secondary structure derived by DSSR')


def p():
    """A shortcut for putting a seq at the bottom. Pretty cool for screenshots with names of objects.

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/p.png
    """
    cmd.set("seq_view_format", 4)
    cmd.set("seq_view", 1)
    cmd.set("seq_view_location", 1)
    cmd.set("seq_view_overlay", 1)


def rna_cartoon():
    """http://www-cryst.bioc.cam.ac.uk/members/zbyszek/figures_pymol

    .. image:: ../pngs/rna_cartoon.png
    """
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_ladder_mode", 1)
    cmd.set("cartoon_nucleic_acid_mode", 4)
    cmd.set("cartoon_ring_transparency", 0.5)


def rp17():
    """Color-coding for secondary structure elements for the RNA Puzzle 17.

    For the variant::

         CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAAAUAUCAGGUGCAA
         ((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))
         # len 62-nt

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/rna.png
    """
    txt = """color forest, resi 1-5+12-16; # p1
 color magenta, resi 6-11+34-39;
 color grey, resi 17-24;
 color marine, resi 25-28+59-62;
 color deepblue, resi 29-33+40-42;
 color orange, resi 44-47+48-56;
 color yellow, resi 57-58;
 color red, resi 19+20+21;
"""
    color_by_text(txt)


def rp172():
    """Color-coding for secondary structure elements for the RNA Puzzle 17.

    For the variant::

         CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAUCAGGUGCAA
         ((((.[[[[[[.))))........((((.....]]]]]]...((((()))))..))))
         # len 58-nt

    See rp17()
    """

    txt = """color forest, resi 1-5+12-16; # p1
 color magenta, resi 6-11+34-39
 color grey, resi 17-24
 color marine, resi 25-28+55-58
 color deepblue, resi 29-33+40-42;
 color orange, resi 43-47+48-52;
 color yellow, resi 53-54;
 color red, resi 19+20+21;
"""
    color_by_text(txt)

def color_aa_types():
    """Color aminoacides types like in Cider (http://pappulab.wustl.edu/CIDER/)"""
    txt = """
color gray70, resn Ala+Ile+Leu+Met+Phe+Trp+Val #hydrophobic
color yellow, resn Tyr+Trp #aromatic
color blue, resn Arg+Lys+His # positive
color forest, resn GLN+SER+GLY+thr
color pink, resn PRO # pro
color red, resn GLU+asp # """
    print("""color (according to) amino-acids types)
hydrohobic (gray)  Ala+Ile+Leu+Met+Phe+Trp+Val
aromatic (yellow) Tyr+Trp
positive (blue)  Arg+Lys+His
polar (forest) Gln+Ser+Glu+Thr
negative (red) Glu+Asp
prolina ;) (pink) Pro""")
    color_by_text(txt)


def color_obj(rainbow=0):

        """
        stolen from :)
AUTHOR
        Gareth Stockwell

USAGE
        color_obj(rainbow=0)

        This function colours each object currently in the PyMOL heirarchy
        with a different colour.  Colours used are either the 22 named
        colours used by PyMOL (in which case the 23rd object, if it exists,
        gets the same colour as the first), or are the colours of the rainbow

        """

        # Process arguments
        rainbow = int(rainbow)

        # Get names of all PyMOL objects
        obj_list = cmd.get_names('objects')

        if rainbow:

           print "\nColouring objects as rainbow\n"

           nobj = len(obj_list)

           # Create colours starting at blue(240) to red(0), using intervals
           # of 240/(nobj-1)
           for j in range(nobj):
              hsv = (240-j*240/(nobj-1), 1, 1)
              # Convert to RGB
              rgb = hsv_to_rgb(hsv)
              # Define the new colour
              cmd.set_color("col" + str(j), rgb)
              print obj_list[j], rgb
              # Colour the object
              cmd.color("col" + str(j), obj_list[j])

        else:
           # List of available colours
           colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
           'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', \
           'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',    \
           'wheat', 'white', 'grey' ]
           ncolours = len(colours)

           # Loop over objects
           i = 0
           for obj in obj_list:
              print "  ", obj, colours[i]
              cmd.color(colours[i], obj)
              i = i+1
              if(i == ncolours):
                 i = 0


def names():
    # Get names of all PyMOL objects
    obj_list = cmd.get_names('objects')
    for o in obj_list:
        print(o)


def color_rbw(rainbow=0):
        """
        similar to color_obj() but this time colors every obect as rainbow
        """
        rainbow = int(rainbow)

        # Get names of all PyMOL objects
        obj_list = cmd.get_names('objects')

        if rainbow:

           print "\nColouring objects as rainbow\n"

           nobj = len(obj_list)

           # Create colours starting at blue(240) to red(0), using intervals
           # of 240/(nobj-1)
           for j in range(nobj):
              hsv = (240-j*240/(nobj-1), 1, 1)
              # Convert to RGB
              rgb = hsv_to_rgb(hsv)
              # Define the new colour
              cmd.set_color("col" + str(j), rgb)
              print obj_list[j], rgb
              # Colour the object
              cmd.color("col" + str(j), obj_list[j])
        else:
           colours = ['rainbow']
           ncolours = len(colours)

           # Loop over objects
           i = 0
           for obj in obj_list:
              print "  ", obj, colours[i]
              cmd.spectrum('count', colours[i], obj)
#              cmd.color(colours[i], obj)
              i = i+1
              if(i == ncolours):
                 i = 0

def ino():
    """Sphare and yellow inorganic, such us Mg.

    .. image:: ../../rna_pdb_tools/utils/PyMOL4RNA/doc/ion.png"""
    cmd.show("spheres", "inorganic")
    cmd.set('sphere_scale', '0.25', '(all)')
    cmd.color("yellow", "inorganic")


def spli():
    AllObj = cmd.get_names("all")
    for name in AllObj:
        if 'Exon' in name or 'exon' in name:
            cmd.color('yellow', name)
        if 'Intron' in name or 'intron' in name:
            cmd.color('magenta', name)
        if name.startswith("U2_snRNA"):
            cmd.color('green', name)
        if name.startswith("U5_snRNA"):
            cmd.color('blue', name)
        if name.startswith("U4_snRNA"):
            cmd.color('orange', name)
        if name.startswith("U4_snRNA"):
            cmd.color('red', name)
    cmd.color("red",'resn rG+G and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2')
    cmd.color("forest",'resn rC+C and name n1+c2+o2+n3+c4+n4+c5+c6')
    cmd.color("orange",'resn rA+A and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2')
    cmd.color("blue",'resn rU+U and name n3+c4+o4+c5+c6+n1+c2+o2')


def _spli():
    """
    # this trick is taken from Rhiju's Das code
    color red,resn rG+G and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2
    color forest,resn rC+C and name n1+c2+o2+n3+c4+n4+c5+c6
    color orange, resn rA+A and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2
    color blue, resn rU+U and name n3+c4+o4+c5+c6+n1+c2+o2

    #
    #cmd.color("yellow", "*intron*")
    #cmd.color("yellow", "*exon*")

    #cmd.show("spheres", "inorganic")
    #cmd.color("yellow", "inorganic")
    """
    cmd.color("orange", "U4_snRNA*")
    cmd.color("red", "U6_snRNA*")
    cmd.color("blue", "U5_snRNA*")
    cmd.color("green", "U2_snRNA*")
    cmd.color("red",'resn rG+G and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2')
    cmd.color("forest",'resn rC+C and name n1+c2+o2+n3+c4+n4+c5+c6')
    cmd.color("orange",'resn rA+A and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2')
    cmd.color("blue",'resn rU+U and name n3+c4+o4+c5+c6+n1+c2+o2')


def rgyration(selection='(all)', quiet=1):
    '''

[PyMOL] RES: radius of gyration
From: Tsjerk Wassenaar <tsjerkw@gm...> - 2011-03-31 14:07:03
https://sourceforge.net/p/pymol/mailman/message/27288491/
DESCRIPTION

    Calculate radius of gyration

USAGE

    rgyrate [ selection ]
 :::warning:::
 if nothing is selected  function is calculating radius of gyration for all pdbs in current Pymol session
    '''
    quiet = int(quiet)
    model = cmd.get_model(selection).atom
    x = [i.coord for i in model]
    mass = [i.get_mass() for i in model]
    xm = [(m*i,m*j,m*k) for (i,j,k),m in izip(x,mass)]
    tmass = sum(mass)
    rr = sum(mi*i+mj*j+mk*k for (i,j,k),(mi,mj,mk) in izip(x,xm))
    mm = sum((sum(i)/tmass)**2 for i in izip(*xm))
    rg = math.sqrt(rr/tmass - mm)
    if not quiet:
        print "Radius of gyration: %.2f" % (rg)
    return rg


# main code #
try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
else:
    print('   PyMOL4RNA (rna-pdb-tools)  ')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('Quickref: ')
    print('  alter (sele), chain="B" ')
    print('  alter (sele), resv -= 4')
    print('  alter (chain B), resv -= 44 ')
    print('set dash_color, red; set dash_width, 4')
    print('p - prepare seq for printing')
    print('rp - rna present, object names only click to get compact legend')
    print('rp17')
    print('get_pdb')
    print('rna_cartoon')
    print('rs')
    print('rcomp')
    print('color_obj')
    print('color_rbw')
    print('aa')
    print("""cspli - color snRNAs of the spliceosome:
    green: U2,  blue: U5, red:U6, orange:U2""")
    print('RNA_PDB_TOOLS env variable used: ' + RNA_PDB_TOOLS)

    cmd.extend('rp17', rp17)
    cmd.extend('rp', rp)
    cmd.extend('p', p)
    cmd.extend('get_pdb', get_pdb)
    cmd.extend('get_seq', get_seq)
    cmd.extend('rna_cartoon', rna_cartoon)
    cmd.extend('rs', rs)
    cmd.extend('ino', ino)
    cmd.extend('rcomp', rcomp)
    cmd.extend('color_obj', color_obj)
    cmd.extend('color_rbw', color_rbw)
    cmd.extend('aa', align_all)
    cmd.extend('ss', ss)
    cmd.extend('ss_all', ss_all)
    cmd.extend('clarna', clarna)
    cmd.extend("rgyration", rgyration)
    cmd.extend("spli", spli)

    cmd.extend('color_aa_types', color_aa_types)

    cmd.extend('names', names)

    # set dash lines
    cmd.set('dash_color', 'red')
    cmd.set('dash_width', 4)

    print('###########################')
    print('PYMOL4RNA loading .... [ok]')
    print('###########################')
