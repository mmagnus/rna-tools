https://pymolwiki.org/index.php/Category:CGO

# Draw a point

    pseudoatom pt1, pos=[x1, y1, z1]
    pseudoatom pt2, pos=[x2, y2, z2]
    distance /pt1, /pt2
    set dash_gap, 0
    
    You can also create CGO lines between two points, if you want.

https://sourceforge.net/p/pymol/mailman/message/25795427/

# Draw a circle

![](docs/circle.png)

`cgoCircle(55,35,41,.1)` on GC.pdb

https://pymolwiki.org/index.php/CgoCircle
