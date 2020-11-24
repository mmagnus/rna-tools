def rp17():
    """Color-coding for secondary structure elements for the RNA Puzzle 17.

    For the variant::

         CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAAAUAUCAGGUGCAA
         ((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))
         # len 62-nt

    .. image:: ../../rna_tools/tools/PyMOL4RNA/doc/rna.png
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

def rp17csrv():
    """Color-coding for secondary structure elements for the RNA Puzzle 17.

    For the variant::

         CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAAAUAUCAGGUGCAA
         ((((.[[[[[[.))))........((((.....]]]]]]...(((((....)))))..))))
         # len 62-nt

    .. image:: ../../rna_tools/utils/PyMOL4RNA/doc/rna.png
    """
    txt = """color forest, resi 1-5+12-16; # p1
 color magenta, resi 6-11+34-39;
 color grey, resi 17-24;
 color marine, resi 25-28+59-62;
 color deepblue, resi 29-33+40-42;
 color orange, resi 44-47+48-56;
 color yellow, resi 57-58;
 color red, resi 5+19+20+21+31+32+33+40+41+42
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


def rp06():
  txt = """color black, all
  color pink, resi 2-10+163-170
  color grey, resi 12-33
  color green, resi 40-41
  color green, resi 161-162
  color orange, resi 45-61
  color green, resi 64-73
  color blue, resi 74-155
  color cyan, resn B1Z"""
  for t in txt.split('\n'):
    color, resi = t.replace('color ', '').split(',')
    print(color, resi)
    cmd.color(color.strip(), resi.strip())


def rp14():
  """color black; # everything
 color blue, resi 1-5+55-59; # p1
 color green, resi 7-11+16-20; # p2
 color magenta, resi 23+60; # pk
 color yellow, resi 29-34+45-50; # p3
 color grey, resi 24-28+51-54; # e-loop
 color red, resi 6+21+22+24+25+28+52+54; # higly conserved"""
 #color blue, resi 5+55


  txt ="""color black, all
 color red, resi 1-5+55-59
 color blue, resi 1-5+55-59; # p1
 color green, resi 7-11+16-20
 color magenta, resi 23+60
 color yellow, resi 29-34+45-50
 color grey, resi 24-28+51-54
 color red, resi 6+21+22"""
  for t in txt.split('\n'):
    color, resi = t.replace('color ', '').split(',')
    print(color, resi)
    cmd.color(color.strip(), resi.strip())

def rp14s():
  """color with Baker's SHAPE data for rp14!"""
  txt = """
   color yellow, resi 12-15+25-29+35-44
   color red, resi 21-24+53+54+60
  """
  color_by_text(txt)

cmd.extend('rp17', rp17)
cmd.extend('rp17csrv', rp17csrv)
