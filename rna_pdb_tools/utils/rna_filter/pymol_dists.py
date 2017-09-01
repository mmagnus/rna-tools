"""show_dists - show distances in PyMOL

.. image:: ../pngs/getdists.png

Usage::

    PyMOL>show_dists([[1,2]])
    1, 2, 3.41

"""


def show_dists(interactions):  # l=([1,2], [3,4])
    for i in interactions:
        a = "////" + str(i[0]) + "/C2"
        b = "////" + str(i[1]) + "/C2"
        x = cmd.distance('d' + str(i[0]) + '-' + str(i[1]), "(" + a + ")", "(" + b + ")")
        print("%i, %i, %s" % (i[0], i[1], str(round(x, 2))))  # why? to have 2.xx
        # mode, 4)


try:
    from pymol import cmd
except:
    print('PyMOL required')
else:
    print('show_dists([[1,2]]) loaded')
    cmd.extend('show_dists', show_dists)
