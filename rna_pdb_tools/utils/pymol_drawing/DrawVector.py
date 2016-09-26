"""https://pymolwiki.org/index.php/CGOCylinder"""

def draw_vector(x1,y1,z1,x2,y2,z2):
    radius = 0.1
    r1,g1,b1 = 0,0,1 # color (blue)
    r2,g2,b2 = 1,0,0 # color (red)    
    cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "vector" )

cmd.extend("draw_vector", draw_vector)    
