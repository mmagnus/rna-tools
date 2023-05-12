#!/usr/bin/env python
"""
Quick reference:

- clarna: show contacts classification of the selected fragment based on ClaRNA
- ss: show secondary structure of the selection based on py3dna.py (3DNA (Lu, Olson 2003))
- ss_all: the same as ss() but for all objects
- pdbsrc: show PDB content (source) of selection.
- seq: show sequence of the selection
- ino: represent ions as sphare and yellow inorganic, such us Mg
- p: shortcut for putting a seq at the bottom. Pretty cool for screenshots with names of objects
- spli: color snRNA of the spliceosome and bases according to identity U(blue), A(orange), G(red), C(forest)
- rp: @todo
- rs: @todo
- rib: @todo
- select br. all within 12 of resi 574

If you want more, read for interesting functions <https://daslab.stanford.edu/site_data/docs_pymol_rhiju.pdf>

Tips; cmd.do

"""
# imports
import tempfile
import math
import subprocess
import os
import sys

import getpass
user = getpass.getuser()

import os
import sys
from rna_tools.tools.PyMOL4RNA.libs.show_contacts import show_contacts
from rna_tools.tools.PyMOL4RNA.libs.get_raw_distances import get_raw_distances


from pymol import cmd
# axes.py
from pymol.cgo import *
from pymol.vfont import plain

try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)

import imp
try:
    from rna_tools.rna_tools_lib import RNAStructure
    from rna_tools.tools.PyMOL4RNA import code_for_spl
    imp.reload(code_for_spl)
    print('code_for_spl loaded...')
    from rna_tools.tools.PyMOL4RNA import code_for_color_spl
    imp.reload(code_for_color_spl)
    print('code_for_color_spl loaded...')
    import rna_tools
    RNA_TOOLS_PATH = rna_tools.rna_tools_lib.get_rna_tools_path()
    sys.path.insert(0, RNA_TOOLS_PATH) # '/Users/magnus/work/src/rna-tools/rna_tools/tools/PyMOL4RNA')#
    #(os.path.abspath(os.path.dirname(__file__))))
    from rna_tools.tools.PyMOL4RNA import PyMOL4Spliceosome
    imp.reload(PyMOL4Spliceosome)
except ImportError:
    print("rna_tools lib is missing")
    RNA_TOOLS_PATH = ''
    
try:
    RNA_TOOLS_PATH
    EXECUTABLE
except NameError:
    EXECUTABLE="/bin/zsh"
    SOURCE=""

try:
    cmd.set('cartoon_gap_cutoff', 0)
except:
    pass

def exe(cmd, verbose=False):
    """Helper function to run cmd. Using in this Python module."""
    if verbose: print('cmd:' + cmd)
    o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         executable=EXECUTABLE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err

def spla():
    cmd.do("color forest, chain A")
    cmd.do("color firebrick, chain B")
cmd.extend('spla', spla)
print('spla - color A and B')

def spln():
    cmd.do("color forest, chain 2")
    cmd.do("color firebrick, chain 6")
cmd.extend('spln', spln)
print('spln - color 2 and 6')
      
def color_protein():
    cmd.do("color blue, resn ARG+LYS+HIS and (sele)")
    cmd.do("color red, resn ASP+GLU and (sele)")
    cmd.do("color green, resn GLY+ALA+VAL+LEU+ILE+MET+PHE and (sele)")
    cmd.do("color yellow, resn TYR+TRP and (sele)")
    cmd.do("color forest, resn SER+THR+CYS+ASN+GLN and (sele)")
    cmd.do("color pink, resn PRO and (sele)")

cmd.extend('cp', color_protein)

def save_transformed(object, file):
    """Saves the molecule with coordinates from the current orientation.

     Args:
        object (string): PyMOL name
        file (string): a file name to output file

    Example::

         PyMOL>save_transformed 6bk8_RNA_only_Oriented, 6bk8_RNA_only_Oriented.pdb

    Source: <https://pymolwiki.org/index.php/Modeling_and_Editing_Structures>
    """
    m = cmd.get_view(0)
    ttt = [m[0], m[1], m[2], 0.0,
           m[3], m[4], m[5], 0.0,
           m[6], m[7], m[8], 0.0,
           0.0,   0.0,  0.0, 1.0]
    cmd.transform_object(object,ttt,transpose=1)
    cmd.save(file,object)


def color_by_text(txt):
    """Helper function used for color-coding based on residue indexes ranges."""
    for t in txt.strip().split('\n'):
        print(t)
        color, resi = t.replace('color ', '').split(',')
        print((color, resi))
        cmd.color(color.strip(), resi.strip())

def cmd_text(txt):
    """Helper function used for color-coding based on residue indexes ranges."""
    for t in txt.strip().split('\n'):
        cmd.do(t)

def delete_all():
    cmd.delete('*')
cmd.extend('dall', delete_all)
print('dall - delete_all')

def rp():
    """Represent your RNA."""
    cmd.hide("sticks", "all")
    cmd.hide("lines", "all")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 2)
    cmd.set("cartoon_ladder_mode", 1)

def show_all_at_once():
    cmd.set('states', 'on')


def grid_on():
    cmd.set('grid_mode', 1)
def grid_off():
    cmd.set('grid_mode', 0)

cmd.extend('gridon', grid_on)
cmd.extend('gridoff', grid_off)
cmd.extend('gn', grid_on)
cmd.extend('gf', grid_off)


def rs():
    """    The function creates super-cool cartoon-like RNA and colors each structure as a rainbow.
    Good to view aligned structures in a grid.

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/rs.png
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
        print("  ", obj, colours[i])
        cmd.spectrum('count', colours[i], obj)
        i = i+1
        if(i == ncolours):
           i = 0


def rx():
    """    The function creates super-cool cartoon-like RNA and colors each structure as a rainbow.
    Good to view aligned structures in a grid.

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/rs.png
    """
    cmd.hide("sticks", "all")
    cmd.hide("lines", "all")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_ring_mode", 0)
    cmd.set("cartoon_ring_finder", 0)
    cmd.set("cartoon_ladder_mode", 0)
    cmd.set("cartoon_ring_transparency", 0.30)

cmd.extend('rx', rx)


def get_intrs_all_vs_all(verbose=True, redundant=True):
    """
    get_intrs_all_vs_all()
    get_raw_distances contacts_all # U6_C-CWC2_C_all

        # if true all vs all (a-b and b-a) then set redundant to True
        # this is sometimes useful if you want to have interactions of b in form of
        # b-a
        # b-c
        # if redundant False then you will have only
        # a-b
        # b-c

    """
    # cmd.delete('contacts')
    objs = cmd.get_names_of_type("object:molecule")
    if verbose:
        print(objs)
    # objs = ['U6_C', 'CWC25_C']
    objs2 = objs.copy()
    for o in objs:
        if not redundant:
            objs2.pop(0)
        if verbose:
            print(' ', objs2)
        for o2 in objs2:
            if o != o2:  # don't compare object to itself
                print(o,'<>',o2)
                results = o + '-' + o2
                if show_contacts(o, o2, results): #, 'contacts) #) # results to U6_C-CWC15_C
                    # if not None
                    p = '/Users/magnus/Desktop/spl-csv/'  # TODO
                    # _all or 
                    get_raw_distances(results + '_all', filename=p + results + '.csv')
            
cmd.extend('get_intrs_all_vs_all', get_intrs_all_vs_all)

def align_all(cycles = 5, filename="_rmsd_.csv"):
    """
    Args:

        cycles (int): maximum number of outlier rejection cycles {default: 5}

    Returns:

        Prints a table of ref vs models with 7 items:

        RaR  RMSD after refinement
        #AA  Number of aligned atoms after refinement
        CoR  Number of refinement cycles
        RbR  RMSD before refinement
        #AbR Number of aligned atoms before refinement
        RS   Raw alignment score
        AR   Number of residues aligned 

        and saves the table to filename as csv

    old version:

          1_solution_0_rpr 1_santalucia_1_rpr 5.60600471496582 958 4 5.763411521911621 974 416.0 46 -- RMSD 5.76  of  46 residues

"""
    molecules = cmd.get_names_of_type("object:molecule")
    ref = molecules.pop(0)
    print("""
    RaR  RMSD after refinement
    #AA  Number of aligned atoms after refinement
    CoR  Number of refinement cycles
    RbR  RMSD before refinement
    #AbR Number of aligned atoms before refinement
    RS   Raw alignment score
    AR   Number of residues aligned 
    """)
    report = []
    header = 'Ref                  Model                RaR  #AA  CoR  RbR  #AbR RS   AR'
    print(header)
    txt = 'Ref,Model,RMSD after refinement,Number of aligned atoms after refinement, Number of refinement cycles, RMSD before refinement, Number of aligned atoms before refinement, Raw alignment score, Number of residues aligned\n'
    for molecule in molecules:
        values = cmd.align(molecule, ref, cycles=cycles)
        l = ([ref[:20].ljust(20), molecule[:20].ljust(20), str(round(values[0], 2)).ljust(4),
            str(round(values[1], 2)).ljust(4),
            str(round(values[2], 2)).ljust(4),
            str(round(values[3], 2)).ljust(4),
            str(round(values[4], 2)).ljust(4),
            str(round(values[5])).ljust(4),
            str(round(values[6], 2)).ljust(4)])
        print(' '.join(l))
        txt += ','.join([x.strip() for x in l]) + '\n'
        report.append([ref, molecule, values[3], values[6]])

    with open(filename, 'w') as f:
        f.write(txt)

cmd.extend('align_all', align_all)


def rmsdx(cycles = 5, matrix_fn = 'matrix.txt'):
    """
    Args:

    cycles (int): refinement cycles of PyMOL align, default: 5
    matrix_fn (string): a file to save the matrix
                        matrix is pretty much saved in space-separated values 
                        so you can load it to pandas with

                        df = pd.read_csv('matrix.txt', sep=' ', index_col=False)
                        df = df.set_index(df.columns)
                        print(df)
                                      Bact_7DCO_S  Bact_5gm6_S
                        Bact_7DCO_S         0.000        0.562

    Returns:

    string: matrix 
            and matrix_fn file ;-)

    """
    models = cmd.get_names_of_type("object:molecule")
    print(' # of models:', len(models))
    
    f = open(matrix_fn, 'w')
    #t = '# ' # for numpy
    t = ''  # for pandas
    for r1 in models:
        # print r1,
        t += str(r1) + ' '  # here ' ' could be changed to , or \t
    t = t.strip() + '\n'

    c = 1
    for r1 in models:
        for r2 in models:
            if r1 == r2:
                rmsd = 0
            else:
                print(r1, r2)
                values = cmd.align(r1, r2, cycles=cycles)
                # RaR [1]       RbR [3]
                # RaR  #AA  CoR               RbR  #AbR RS   AR'
                # (0.668652355670929, 241, 5, 1.1646124124526978, 293, 199.0, 38)
                rmsd = round(values[0], 3)
            t += str(rmsd) + ' '
        #print('...', c, r1)
        c += 1
        t += '\n'

    f.write(t)
    f.close()

    print(t.strip())  # matrix
    return t

cmd.extend('rmsdx', rmsdx)

def save_all(dir=''):
    """save_all molecule objects as pdb files. Use `cd` to get to the right folder
       or use dir argument"""
    if dir:
        dir += '/'
    molecules = cmd.get_names_of_type("object:molecule")
    for molecule in molecules:
        print('Saving %s ...' % molecule)
        cmd.save(dir + molecule + '.pdb', molecule)

cmd.extend('save_all', save_all)

def save_each_object(folder='', prefix=''):
    """

    Usage::

        save_each_object /Users/magnus/work/spliceosome/PyMOL4Spliceosome/chains/yB_5zwo, yB_5zwo_

        p = 'yP_6exn' # yP_5ylz' #yI_5y88' # yE_6n7r'
        pth = '/Users/magnus/work/spliceosome/PyMOL4Spliceosome/chains/'
        save_each_object(pth + p, p + '_')

    See the application here <https://github.com/mmagnus/PyMOL4Spliceosome/releases/tag/v0.32>

    .. todo:: add some way to select which objects to use
    """
    obj_list = cmd.get_names('objects')
    for o in obj_list:
        if folder:
            folder += '/'
        fn = folder + prefix.strip() + o.strip() + '.pdb'
        cmd.save(fn, o)

cmd.extend('save_each_object', save_each_object)

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
        print("  ", obj, colours[i])
        cmd.color(colours[i], obj)
        i = i+1
        if(i == ncolours):
           i = 0

def da():  # donors accepters
    t = """
    h_add;
    set sphere_scale, 0.2, (all)
    set sphere_transparency, 0
    #color blue, donors;
    #color green, acceptors;
    show sphere, donors;
    show sphere, acceptors;
    color gray,  name H*; # color atoms from white to gray
    """
    cmd.do(t)
cmd.extend('da', da)
    
def pdb():
    """Get PDB content of selection.

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/pdb.png"""
    tmpfn = '/tmp/pymol_get_pdb.pdb'
    cmd.save(tmpfn, '(sele)')
    s = RNAStructure(tmpfn)
    for l in s.lines:
        print(l)


def x3dna():
    f = tempfile.NamedTemporaryFile(delete=False) # True)
    #cmd.save(f.name + '.pdb', '(backbone_)')
    cmd.save(f.name + '.pdb', '(sele)')
    out, err = exe("rna_x3dna.py --show-log " + f.name + ".pdb ")
    print('\n'.join(out.split('\n')[1:]))  # to remove first line of py3dna /tmp/xxx
    if err:
        print(err)
    f.close()


def clarna2():
    """Get contacts classification of the selected fragment based on ClaRNA (for each object).

    .. image:: ../../rna_tools/tools/PyMOL4RNA/doc/clarna.png
    """
    objs = cmd.get_names("objects")
    for name in objs[:]:
        print(name + ' ' + '-' * (70 - len(name)))
        f = tempfile.NamedTemporaryFile(delete=False) # True)
        #cmd.save(f.name + '.pdb', '(backbone_)')
        cmd.save(f.name + '.pdb', '(sele) and "' + name + '"')
        CLARNA_RUN = 'rna_clarna_run.py'
        #cmdline = #SOURCE + " && " +
        cmdline = CLARNA_RUN + " -ipdb " + f.name + '.pdb -bp+stack'
        print(cmdline)
        out, err = exe(cmdline)
        print('\n'.join(out.split('\n')[1:]))  # to remove first line of py3dna /tmp/xxx
        if err:
            print(err)
        f.close()


def clarna(selection, folder:str=''):
    """Get contacts classification of the selected fragment based on ClaRNA (for each object).

    Args:

      folder (str): The path to save temporary files, by default they are save to system tmp

    Example::

        PyMOL>clarna sele
        rna_clarna_run.py -ipdb /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmp1h_bwvtx.pdb -bp+stack
        chains:  X 15 16
        X   15   X   16          bp A A                      ><   0.9427

    .. image:: ../../rna_tools/tools/PyMOL4RNA/doc/clarna.png
    """
    f = tempfile.NamedTemporaryFile(delete=False) # True)
    if not folder:
        output = f.name + '_clarna.pdb'
    else:
        output= folder + os.sep + os.path.basename(f.name) + '_clarna.pdb'

    cmd.save(output, selection)
    CLARNA_RUN = 'rna_clarna_run.py'
    cmdline = CLARNA_RUN + " -ipdb " + output + ' -bp+stack'
    print(cmdline)
    out, err = exe(cmdline)
    print('\n'.join(out.split('\n')[1:]))  # to remove first line of py3dna /tmp/xxx
    if err:
        print(err)
    cmd.load(output)
    f.close()


def seq(selection):
    """Get sequence of the selected fragment using ``rna_pdb_tools.py --get_seq ``.

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/ss.png
    """
    if selection.strip() == "*":
        AllObj = cmd.get_names("all")
        # print AllObj
        for name in AllObj[:]:
            if not name.startswith('_align'):
                f = tempfile.NamedTemporaryFile(delete=False) # True)
                f.name = f.name + '.pdb'
                cmd.save(f.name, name)
                cmdline = 'rna_pdb_tools.py --color-seq --get-seq ' + f.name
                out, err = exe(cmdline)
                if out:
                    print('> ' + name)
                    print('\n'.join(out.split('\n')[2:]))  # to remove first line of py3dna /tmp/xxx
                    # hide this line: is >tmpGCszi7 nts=4 [tmpGCszi7] -- secondary structure derived by DSSR
                if err:
                    # print(err)
                    pass
                f.close()
    else:
        f = tempfile.NamedTemporaryFile(delete=False)
        selection = strip_selection_name(selection)
        input = os.path.dirname(f.name) + os.sep +  selection + '.pdb'
        cmd.save(input, selection)
        cmdline = 'rna_pdb_tools.py  --color-seq --get-seq ' + input
        # print(cmdline)
        out, err = exe(cmdline)
        print(out)
        if err:
            print(err)
        f.close()

def seqsel():
    """Get sequence of the selected fragment using ``rna_pdb_tools.py --get_seq ``.
    """
    f = tempfile.NamedTemporaryFile(delete=False)
    selection = '(sele)'
    input = os.path.dirname(f.name) + os.sep +  '_sele.pdb'
    cmd.save(input, selection)

    cmdline = 'rna_pdb_tools.py --get-seq ' + input
    print(cmdline)
    out, err = exe(cmdline)
    print(out)
    if err:
        print(err)
    f.close()


def ss(selection):
    """Get Secondary Structure of (sele) based on py3dna.py.

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/ss.png
    """
    f = tempfile.NamedTemporaryFile(delete=False)
    output = os.path.dirname(f.name) + os.sep +  selection + '.pdb'
    cmd.save(output, '(sele)')

    cmdline = 'rna_x3dna.py ' + output
    out, err = exe(cmdline)
    print('\n'.join(out.split('\n')[2:]))  # to remove first line of py3dna /tmp/xxx
    if err:
        print(err)
    f.close()



def rtrun(cmd, selection, suffix):
    f = tempfile.NamedTemporaryFile(delete=False) # True)
    output = os.path.dirname(f.name) + os.sep +  selection + '.pdb'
    output2 = os.path.dirname(f.name) + os.sep + selection + '_mut.pdb'
    exe(cmdline)
    print(cmdline)
    cmd.save(output, selection)

    # 'A:1A+2A+3A+4A' 

def mutate(mutation, selection):
    """
    """
    f = tempfile.NamedTemporaryFile(delete=False) # True)
    output = os.path.dirname(f.name) + os.sep +  selection + '.pdb'
    output2 = os.path.dirname(f.name) + os.sep + selection + '_mut.pdb'
    cmdline = "rna_pdb_tools.py --mutate " + mutation + ' ' + output + ' > ' + output2
    print(cmdline)
    exe(cmdline)
    cmd.load(output2)

cmd.extend('mutate', mutate)


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
            out, err = exe(RNA_TOOLS_PATH + '/bin/rna_x3dna.py ' + f.name)
            print('\n'.join(out.split('\n')[2:]))  # to remove first line of py3dna /tmp/xxx
            # hide this line: is >tmpGCszi7 nts=4 [tmpGCszi7] -- secondary structure derived by DSSR
            if err:
                print(err)
            f.close()
    print('-- secondary structure derived by DSSR')


def p():
    """A shortcut for putting a seq at the bottom. Pretty cool for screenshots with names of objects.

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/p.png
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

           print("\nColouring objects as rainbow\n")

           nobj = len(obj_list)

           # Create colours starting at blue(240) to red(0), using intervals
           # of 240/(nobj-1)
           for j in range(nobj):
              hsv = (240-j*240/(nobj-1), 1, 1)
              # Convert to RGB
              rgb = hsv_to_rgb(hsv)
              # Define the new colour
              cmd.set_color("col" + str(j), rgb)
              print(obj_list[j], rgb)
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
              print("  ", obj, colours[i])
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

           print("\nColouring objects as rainbow\n")

           nobj = len(obj_list)

           # Create colours starting at blue(240) to red(0), using intervals
           # of 240/(nobj-1)
           for j in range(nobj):
              hsv = (240-j*240/(nobj-1), 1, 1)
              # Convert to RGB
              rgb = hsv_to_rgb(hsv)
              # Define the new colour
              cmd.set_color("col" + str(j), rgb)
              print(obj_list[j], rgb)
              # Colour the object
              cmd.color("col" + str(j), obj_list[j])
        else:
           colours = ['rainbow']
           ncolours = len(colours)

           # Loop over objects
           i = 0
           for obj in obj_list:
              print("  ", obj, colours[i])
              cmd.spectrum('count', colours[i], obj)
#              cmd.color(colours[i], obj)
              i = i+1
              if(i == ncolours):
                 i = 0


def strip_selection_name(selection_name):
    """Quick function: (sele) -> sele"""
    return selection_name.replace('(', '').replace(')', '')


def edges(selection):
    """Save selection into a file in a temp folder and run rna_draw_edges.py on it and load it into this session"""
    my_view = cmd.get_view()
    f = tempfile.TemporaryDirectory()
    tmpf = f.name + os.sep + strip_selection_name(selection) + '.pdb'
    outf = f.name + '/output.py'
    cmd.save(tmpf, selection)
    cmdline = '/Users/magnus/miniconda3/bin/rna_draw_edges.py --name %s %s > %s' % (strip_selection_name(selection), tmpf, outf)
    print(cmdline)
    out, err = exe(cmdline)
    if err:
        print(err)
    cmd.load(outf)
    cmd.set_view(my_view)

cmd.extend('edges', edges)


def ino():
    """Sphare and yellow inorganic, such us Mg.

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/ion.png"""
    cmd.show("spheres", "inorganic")
    cmd.set('sphere_scale', '0.25', '(all)')
    #cmd.set('sphere_scale', '1', '(all)')
    cmd.color("yellow", "inorganic")

mapping = [[u'PRP8', 'A', u'skyblue'], [u'BRR2', 'B', u'grey60'], [u'BUD31', 'C', u'dirtyviolet'], [u'CEF1', 'D', u'raspberry'], [u'CLF1', 'E', u'raspberry'], [u'CWC15', 'F', u'dirtyviolet'], [u'CWC16/YJU2', 'G', u'lightteal'], [u'CWC2', 'H', u'ruby'], [u'CWC21', 'I', u'violetpurple'], [u'CWC22', 'J', u'bluewhite'], [u'CWC25', 'K', u'deepteal'], [u'Intron', 'L', u'black'], [u'ISY1', 'M', u'dirtyviolet'], [u'LEA1', 'N', u'palegreen'], [u'Msl1', 'O', u'palegreen'], [u'PRP45', 'P', u'lightpink'], [u'PRP16', 'Q', u'smudge'], [u'CDC40\xa0(PRP17, SLU4, XRS2)', 'R', u'dirtyviolet'], [u'PRP19 (PSO4)', 'S', u'grey70'], [u'PRP46', 'T', u'lightblue'], [u'SLT11/ECM2', 'U', u'chocolate'], [u'SNT309', 'V', u'grey70'], [u'SNU114', 'W', u'slate'], [u'SYF2', 'X', u'brightorange'], [u'SYF1', 'Y', u'brightorange'], [u'U2', 'Z', u'forest'], [u'U5', 'a', u'density'], [u'U5_SmRNP', 'b', u'deepblue'], [u'U6', 'c', u'firebrick'], [u'Intron', 'r', u'grey50'], [u'Exon', 'z', u'yellow'], [u'exon-3', 'y', u'yellow'], [u'exon-5', 'z', u'yellow'], [u'PRP4 ', 'd', u'grey50'], [u'PRP31', 'e', u'grey50'], [u'PRP6', 'f', u'grey50'], [u'PRP3', 'g', u'grey50'], [u'DIB1', 'h', u'grey50'], [u'SNU13', 'i', u'grey50'], [u'LSM8', 'j', u'grey50'], [u'LSM2', 'k', u'grey50'], [u'LSM3', 'l', u'grey50'], [u'LSM6', 'm', u'grey50'], [u'LSM5', 'n', u'grey50'], [u'LSM7', 'o', u'grey50'], [u'LSM4', 'p', u'grey50'], [u'SNU66', 'q', u'grey50'], [u'RNA (intron or U6 snRNA)', 'r', u'grey50'], [u'5EXON', 's', u'grey50'], [u'BUD13', 't', u'grey60'], [u'CLF2', 'u', u'rasberry'], [u'Cus1', 'v', u'palegreen'], [u'CWC24', 'w', u'grey60'], [u'CWC27', 'x', u'grey60'], [u'HSH155', '1', u'smudge'], [u'HSH49', '2', u'sand'], [u'PML1', '3', u'grey60'], [u'PRP11', '4', u'palegreen'], [u'PRP2', '5', u'palegreen'], [u'RDS3', '6', u'palegreen'], [u'RSE1', '7', u'smudge'], [u'SNU17', '8', u'grey60'], [u'Ysf3', '9', u'palegreen'], [u'cwc23', 'd', u'grey50'], [u'SPP382\xa0(CCF8, NTR1)', 'e', u'grey50'], [u'NTR2', 'f', u'grey50'], [u'PRP43', 'g', u'grey50'], [u'SMB1', 'h', u'grey50'], [u'SME1', 'i', u'grey50'], [u'SMX3', 'j', u'grey50'], [u'SMX2\xa0(SNP2)', 'k', u'grey50'], [u'SMD3', 'l', u'grey50'], [u'SMD1', 'm', u'grey50'], [u'SMD2', 'n', u'grey50'], [u'PRP22', 'o', u'grey50'], [u'PRP18', 'p', u'grey50'], [u'SLU7', 'q', u'grey50'], [u'SMF', 'd', u'grey50'], [u'SMG', 'e', u'grey50'], [u'PRP9', 'f', u'grey50'], [u'PRP21', 'g', u'grey50'], [u'SNU23', 'r', u'grey50'], [u'PRP38', 's', u'grey50'], [u'SPP381', 'w', u'grey50']]


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
    from itertools import izip
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
        print("Radius of gyration: %.2f" % (rg))
    return rg



def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


def qrnass():
    cmd.save('sele.pdb', '(sele)')
    mini('sele.pdb')


def qrnas():
    subset = "*"
    AllObj=cmd.get_names("all")
    #print AllObj
    for x in AllObj[:]:
      print(x, 'qrnas...')
      #print(AllObj[0],x)
      f = tempfile.NamedTemporaryFile(delete=True)
      #print f.name
      #f.write(XX)
      cmd.save(f.name, x)
      #p = Process(target=mini)
      #p.start()
      mini()
      #cmd.load('out.pdb', 'ref')
      #p.join()
      #print x
      #print '\n'.join(out.split('\n')[1:]) # to remove first line of py3dna /tmp/xxx
      f.close()
      break
    align_all()
    rr()
    cmd.set('grid_mode', 1)


def inspect(name, dont_green=False):
    f = tempfile.NamedTemporaryFile(delete=False)
    cmd.save(f.name + '.pdb', name)
    out, err = exe('rna_pdb_tools.py --inspect ' + f.name + '.pdb')
    if not dont_green:
        cmd.color('green', name)
    for l in out.split('\n'):
         hit = re.findall('chain: (?P<chain>\w+).*\# (?P<residue>\d+)', l)
         if hit: # [('X', '1')]
             cmd.color('red', name + ' and chain ' + hit[0][0] + ' and resi ' + hit[0][1])
    print(out)

cmd.extend('inspect', inspect)
    
def ll():
    cmd.do('hide all')
    cmd.do('show lines')
cmd.extend('ll', ll)
print('ll - show lines only')

def rpr(selection):
    """rpr"""
    f = tempfile.NamedTemporaryFile(delete=False)
    input = os.path.dirname(f.name) + os.sep +  selection + '.pdb'
    cmd.save(input, selection)

    output = os.path.dirname(f.name) + os.sep +  selection + '_rpr.pdb'

    out, err = exe('rna_pdb_tools.py --rpr ' + input + ' > ' + output)
    cmd.load(output)
    print(out)
cmd.extend('rpr', rpr)

def diff(selection, selection2):
    """rpr"""
    f = tempfile.NamedTemporaryFile(delete=False)
    input = os.path.dirname(f.name) + os.sep +  selection + '.pdb'
    cmd.save(input, selection)
    output = os.path.dirname(f.name) + os.sep +  selection2 + '.pdb'
    cmd.save(input, selection)
    cmdline = 'diffpdb.py ' + input + ' ' + output + ' &'
    print(cmdline)
    #os.system(cmdline)
    exe(cmdline)

cmd.extend('diff', diff)

def mini(f):  # min with qrna

    #os.system('/home/magnus/opt/qrnas/QRNA02/QRNA -i ' + f + ' -c /home/magnus/opt/qrnas/QRNA02/configfile.txt -o out.pdb')
    os.system('~/opt/qrnas/QRNA02/QRNA -i ' + f + ' -c ~/opt/qrnas/QRNA02/configfile.txt -o out.pdb')
    cmd.delete('mini')
    cmd.load('out.pdb', 'mini')
    print('end')


def reload():
    """Reload ~/.pymolrc and all included there packages (e.g. with run <foo.py>)"""
    cmd.do('@ ~/.pymolrc')


def rlabel():
    cmd = "n. C1'", '"%s %s" % (resn, resi)'
    print('label ' + cmd)
    cmd.label(cmd)



def sav(name): #sav
    # cmd.bg_color( "white" )
    tf = tempfile.NamedTemporaryFile(delete=False)
    fn = tf.name + '.png'
    tf = tempfile.NamedTemporaryFile(delete=False)
    cfn = tf.name + '.png'

    psefn = '~/Desktop/' + name + '.pse'
    cmd.save(psefn)

    cmd.save(fn)

    cmdline= "/opt/homebrew/bin/convert " + fn + " -gravity center -crop 3:3 +repage " + cfn
    print(cmdline)
    os.system(cmdline)
    cmdline = '/opt/homebrew/bin/fileicon set ' + psefn + ' ' + cfn
    print(cmdline)
    os.system(cmdline)

    #cmd.png(coverfn, 576,576)
    #cmd.ray(576,576)
cmd.extend('sav', sav)

def hide_rna():
    cmd.hide('(polymer.nucleic)')
cmd.extend('rna-hide', hide_rna)

def show_rna():
    cmd.show('(polymer.nucleic)')
cmd.extend('rna-show', show_rna)

def clr():
    cmd.do('delete *')
cmd.extend('clr', clr)
print('clr - delete all')

def bw():
  """clr - make white bg and structure black"""
  cmd.bg_color( "white" )
  color_by_text('color black, all')
cmd.extend('bw', bw)
print('bw - white bg, black all')

def select_rna():
    cmd.select('polymer.nucleic')
cmd.extend('select-rna', select_rna)

def hide_protein():
    cmd.hide('(polymer.protein)')
#cmd.extend('protein-hide', hide_protein)
#cmd.extend('rp-hide', hide_protein)
def select_protein():
    cmd.select('polymer.protein')
cmd.extend('protein-select', select_protein)

def tmp():
    cmd.save('/home/' + user + '/Desktop/' + tmp + '.png')
    cmd.save('/home/' + user + '/Desktop/' + tmp + '.pse')

def tp(): #tp temp pse
    """tp here"""
    import datetime

    # cmd.bg_color( "white" )
    tf = tempfile.NamedTemporaryFile(delete=False)
    fn = tf.name + '.png'
    tf = tempfile.NamedTemporaryFile(delete=False)
    cfn = tf.name + '.png'

    date = datetime.datetime.today().strftime('%Y-%m-%d.%H%M%S')
    psefn = '~/Desktop/' + date + '.pse'
    cmd.save(psefn)

    cmd.save(fn)

    cmdline= "/opt/homebrew/bin/convert " + fn + " -gravity center -crop 3:3 +repage " + cfn
    print(cmdline)
    os.system(cmdline)
    cmdline = '/opt/homebrew/bin/fileicon set ' + psefn + ' ' + cfn
    print(cmdline)
    os.system(cmdline)

cmd.extend('tp', tp)

def sav_tmp():
    from shutil import copyfile
    import datetime
    try:
        TMP_FOLDER + ' '
    except:
        print("Error: Set up TMP_FOLDER in your ~/.pymolrc, e.g. TMP_FOLDER = '/home/magnus/Desktop/PyMOL/'")
        return

    try:
        os.mkdir(TMP_FOLDER)
    except:
        pass

    date = datetime.datetime.today().strftime('%Y-%m-%d.%H%M%S')
    try:
        fn = TMP_FOLDER +  os.sep + id + '_' + date + '.pse'
    except TypeError:
        fn = TMP_FOLDER +  os.sep + '_' + date + '.pse'
    cmd.save(fn)
    print('Save...' + fn)
    cmd.save(fn.replace('.pse', '.png'))
    copyfile(fn, TMP_FOLDER + '/last.pse')

def load_tmp():
    print('Load...')
    cmd.load(TMP_FOLDER + '/last.pse')


def trim():
    cmd.do("remove solvent")
    cmd.do("remove resn NA")
    
cmd.extend('trim', trim)
    
def get_resi():
    """
    PyMOL>get_resi()
    358+376+288+290+359+383+386+382+289+287+384+357+385+360+377+361
    """
    stored.residues = set()
    cmd.iterate('(sele)', "stored.residues.add(resi)") # ('289') # residue only
    print('+'.join(stored.residues))#[:-1])

    #cmd.iterate('(sele)', "stored.residues.add((chain, resi))") # ('A', '289')
    #print(r, end='+') # ('A', '289')
    #selection = object + ' and index ' + str(index)
    #cmd.iterate('(sele)', 'l.append([resn, resi])')
    #rint(l)

def findN(r):
    c = 'select br. all within ' + str(r) + ' of (sele)'
    cmd.do(c)

def white():
    cmd.set('dash_color', 'black')
    cmd.set('dash_width', 2)
    cmd.bg_color( "white" )
cmd.extend('w', white)

def desc(t='', width=80):
    print()
    print()
    print()
    print(t.center(int(width)))
    print()
    print()
    print()

def s(): # quick save selected to tmp.pdb
    cmd.do('save tmp.pdb, (sele)')
cmd.extend('s', s)
print('s - quick save selected to tmp.pdb')

def draw():
    #cmd.select("name C1'")
    t = """


    set sphere_scale, 0.4, name C1'
    show sphere, name C1'
    color black, name C1'

    remove name OP2+OP1+O5'+P+C5'+O2'+O3'+C4'+C3'+C2'+O4'

    h_add

    select (bound_to name C1') and name H*
    remove (sele)

    set spec_reflect, 0
    # https://pymolwiki.org/index.php/Spec_reflect

    #set sphere_transparency, 0
    #color blue, donors;
    #color green, acceptors;
    #show sphere, donors;
    #show sphere, acceptors;

    color gray,  name H*; # color atoms from white to gray


set dash_color, black
set dash_width, 1; set dash_gap,  0.2
bg white;
#color black;
"""
    cmd_text(t)

cmd.extend('draw', draw)
cmd.extend('dr', draw)

def se(): # save
    cmd.do('save tmp.pdb, (enabled)')
cmd.extend('se', se)
print('se - quick save enabled to tmp.pdb')
      

def axes_big():
    """
    https://pymolwiki.org/index.php/Axes
    """
    cmd.delete('axes')

    # create the axes object, draw axes with cylinders coloured red, green,
    #blue for X, Y and Z

    obj = [
       CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
       CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,
       CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,
       ]

    # add labels to axes object (requires pymol version 0.8 or greater, I
    # believe

    cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,0.,10.],'Z',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])

    # then we load it into PyMOL
    cmd.load_cgo(obj,'axes')

def axes(gradient=False):
    """
    https://pymolwiki.org/index.php/Axes
    """
    cmd.delete('axes')

    # create the axes object, draw axes with cylinders coloured red, green,
    #blue for X, Y and Z

    l = 1
    width = 0.01 #  0.2
    if gradient: # with gradient
        obj = [
           CYLINDER, 0., 0., 0., l, 0., 0., width, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
           CYLINDER, 0., 0., 0., 0., l, 0., width, 1.0, 1.0, 1.0, 0., 1.0, 0.,
           CYLINDER, 0., 0., 0., 0., 0., l, width, 1.0, 1.0, 1.0, 0., 0.0, 1.0,
           ]
    else:
        obj = [
           CYLINDER, 0., 0., 0., l, 0., 0., width, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0., 0., 0., 0., l, 0., width, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0., 0., 0., 0., 0., l, width, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           ]
        
    # add labels to axes object (requires pymol version 0.8 or greater, I
    # believe
    l
    #cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    width = 0.005 # 0.20
    size = 0.05
    cyl_text(obj,plain,[l,0.,0.],'X',width,axes=[[size,0,0],[0,size,0],[0,0,size]])
    cyl_text(obj,plain,[0.,l,0.],'Y',width,axes=[[size,0,0],[0,size,0],[0,0,size]])
    cyl_text(obj,plain,[0.,0.,l],'Z',width,axes=[[size,0,0],[0,size,0],[0,0,size]])

    # then we load it into PyMOL
    cmd.load_cgo(obj,'axes')


cmd.extend('axes', axes)


def axes2(length=0.75):
    """Draw XYZ, no lables here, see axes()
    https://pymolwiki.org/index.php/Axes

    Args:
        lenght: length of axes
    """
    cmd.delete('axes')
    print('Draw axis red, green blue for X, Y and Z')
    w = 0.06 # cylinder width 
    l = float(length) # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter

    obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
           CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
           CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    cmd.load_cgo(obj, 'axes')
cmd.extend('axes2', axes2)

def v(x, y, z, name='v'):
    cmd.delete(name)
    w = 0.01 # cylinder width
    length=0.75
    l = float(length) # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter
    r1,g1,b1 = 1,1,1
    r2,g2,b2 = r1,g1,b1
    obj = [CYLINDER, 0.0, 0.0, 0.0,   x, y, z, w, r1, g1, b1, r2, g2, b2]
    cmd.load_cgo(obj, name)
cmd.extend('v', v)

def ha():
    """
        cmd.do('h_add')
    """
    cmd.do('h_add')
cmd.extend('ha', ha)
    
def hb(): # hydrogen bonds
    cmd.do('contacts *,*')
cmd.extend('hb', hb)


def pn(): # pn
    cmd_text('save ~/Desktop/tmp.png')
cmd.extend('pn', pn) # pn
def quickref():
    print('   PyMOL4RNA (rna-tools)  ')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('Quickref `qr`: ')
    print('  alter (sele), chain="B" ')
    print('  alter (sele), resv -= 4')
    print('  alter (chain B), resv -= 44 ')
    print("  select br. all within 15 of (sele)")
    print("  select br. all within 15 of resi 574")
    print("  select br. all within 15 of resi 377 # O. ihheyensis")
    print('  select br. all within 15 of U6_snRNA and resi 80')
    print('  set dash_color, red; set dash_width, 4')
    print('  p - prepare seq for printing')
    print('  rp - rna present, object names only click to get compact legend')
    print('  rp17')
    print('  rna_cartoon')
    print('  rs')
    print('  rcomp')
    print('  color_obj')
    print('  color_rbw')
    print('  aa')
    print('  findN')
    print('  savt - save_transformed <object>, <file>')
    print(' select br. all within 20 of (sele) #within with aa')
    #print('  spl - color snRNAs of the spliceosome:'
    #    green: U2,  blue: U5, red:U6, orange:U2""")
    print('\_ RNA_TOOLS_PATH env variable used: ' + RNA_TOOLS_PATH)

try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
else:
    quickref()
    #cmd.set_key('CTRL-S', cmd.save, ['/home/magnus/Desktop/tmp.pse'])
    cmd.set_key('CTRL-S', sav_tmp)
    cmd.set_key('CTRL-Z', load_tmp)  # ostatni wrzucam tutaj
    #cmd.load, ['/home/magnus/Desktop/tmp.pse'])
    # main code #

    cmd.extend('quickref', quickref)
    cmd.extend('qr', quickref)

    cmd.extend('rp', rp)
    cmd.extend('p', p)
    cmd.extend('pdb', pdb)
    cmd.extend('seq', seq)
    cmd.extend('seqsel', seqsel)
    cmd.extend('rseq', seq)
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
    cmd.extend('x3dna', x3dna)
    cmd.extend("rgyration", rgyration)
    cmd.extend('rlabel', 'rlabel')

    cmd.extend('reload', reload)
    cmd.extend('rl', reload)

    cmd.extend('color_aa_types', color_aa_types)

    cmd.extend('names', names)

    cmd.extend('findX', findN)

    # set dash lines #hbonds #hydrogen
    cmd.set('dash_color', 'white')
    cmd.set('dash_width', 2)
    cmd.set('cartoon_tube_radius', 0.5)
    
    cmd.extend('save_transformed', save_transformed)
    cmd.extend('savt', save_transformed)
    cmd.extend('show_all_at_once', show_all_at_once)

    cmd.set('ignore_case', 'off')
    #cmd.set('cartoon_ring_mode', '3')
    #cmd.set('cartoon_ring_finder', '2')
    #cmd.extend('spl_select', spl_select)
    #    print('ignore_case made off')
    print('\_ PYMOL4RNA loading .... [ok]')

    cmd.extend('desc', desc)
    #cmd.do('set overlay, 1')

    #### change do desktop
    user = getpass.getuser()
    cw = os.path.abspath(os.getcwd())
    print(cw)
    if cw == '/Users/magnus':
        print('change to Desktop')
        os.chdir('/Users/magnus/Desktop/')
    cw = os.path.abspath(os.getcwd())
    print(cw)
