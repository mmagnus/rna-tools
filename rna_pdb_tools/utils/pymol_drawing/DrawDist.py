"""https://sourceforge.net/p/pymol/mailman/message/25795427/"""

def draw_dist(x1,y1,z1,x2,y2,z2):
    """
    draw_dist(54.729, 28.9375, 41.421, 55.342, 35.3605, 42.745)
    """
    cmd.pseudoatom('pt1', pos=[x1, y1, z1])
    cmd.pseudoatom('pt2', pos=[x2, y2, z2])
    cmd.distance('pt1-pt2', 'pt1','pt2')

cmd.extend("draw_dist", draw_dist)  
