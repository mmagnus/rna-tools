#!/usr/bin/env python

import tempfile
from pymol import cmd
from itertools import izip
import math


from rna_pdb_tools.rna_pdb_tools.rna_pdb_tools import RNAStructure



def color_by_text(txt):
    for t in txt.strip().split('\n'):
        color, resi = t.replace('color ', '').split(',')
        print((color, resi))
        cmd.color(color.strip(), resi.strip())


def rp():
    """RNA like in papers ;-)"""
    cmd.hide("sticks", "all")
    cmd.hide("lines", "all")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 2)
    cmd.set("cartoon_ladder_mode", 1)

    
def rs():
    """RNA like in papers ;-)  - even better :D 
    
    Creates supercool cartoon-like RNA and colors each (and every) structure as a rainbow.
    Good to view aligned structures in a grid.
    
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
    """ """
    tmpfn = '/tmp/pymol_get_pdb.pdb'
    cmd.save(tmpfn, '(sele)')
    s = RNAStructure(tmpfn)
    for l in s.lines:
        print(l)


def __off__ss():
    subset = "*"
    AllObj = cmd.get_names("all")
    # print AllObj
    for x in AllObj[:]:
        # print(AllObj[0],x)
        f = tempfile.NamedTemporaryFile(delete=True)
        # print f.name
        # f.write(XX)
        cmd.save(f.name, x)
        out = subprocess.getoutput('py3dna.py ' + f.name)
        print(x)
        print('\n'.join(out.split('\n')[1:]))  # to remove first line of py3dna /tmp/xxx
        f.close()


def p():
    cmd.set("seq_view_format", 4)
    cmd.set("seq_view", 1)
    cmd.set("seq_view_location", 1)
    cmd.set("seq_view_overlay", 1)


def rna_cartoon():
    """http://www-cryst.bioc.cam.ac.uk/members/zbyszek/figures_pymol

    .. image :: ../pngs/rna_cartoon.png
    """
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_ladder_mode", 1)
    cmd.set("cartoon_nucleic_acid_mode", 4)
    cmd.set("cartoon_ring_transparency", 0.5)


def rp17():
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
def rgyration(selection='(all)', quiet=1):
    '''
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

    cmd.extend('rp17', rp17)
    cmd.extend('rp', rp)
    cmd.extend('p', p)
    cmd.extend('get_pdb', get_pdb)
    cmd.extend('rna_cartoon', rna_cartoon)
    cmd.extend('rs', rs)
    cmd.extend('rcomp', rcomp)
    cmd.extend('color_obj', color_obj)
    cmd.extend('color_rbw', color_rbw)
    cmd.extend('aa', align_all)
    cmd.extend("rgyration", rgyration)
    
    # set dash lines
    cmd.set('dash_color', 'red')
    cmd.set('dash_width', 4)
