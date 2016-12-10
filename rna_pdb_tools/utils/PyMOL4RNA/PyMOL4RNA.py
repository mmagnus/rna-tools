def color_by_text(txt):
  for t in txt.strip().split('\n'):
    color, resi = t.replace('color ', '').split(',')
    print color, resi
    cmd.color(color.strip(), resi.strip())
    
def rp():
    """rna like in papers ;-)"""
    cmd.hide("sticks", "all")
    cmd.hide("lines", "all")
    cmd.show("cartoon", "all")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 2)
    cmd.set("cartoon_ladder_mode", 1)
    
def rp17():
    txt ="""color forest, resi 1-5+12-16; # p1
 color magenta, resi 6-11+34-39
 color grey, resi 17-24
 color marine, resi 25-28+55-58 
 color deepblue, resi 29-33+40-42
 color orange, resi 42-47+48-52;
 color yellow, resi 53-54;
"""
    color_by_text(txt)

cmd.extend('rp17', rp17)
cmd.extend('rp', rp)   
