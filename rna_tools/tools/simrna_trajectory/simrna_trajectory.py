#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SimRNATrajectory module.

SimRNATrajectory / Frame / Residue / Atom
"""

from __future__ import print_function
from collections import deque
import numpy as np
import gc


class SimRNATrajectory:
    """SimRNATrajectory"""

    def __init__(self):
        """Crate SimRNATrajectory object, empty.
        Use
        - load_from_file
        - load_from_string
        to load the data.
        """
        self.frames = []

    def load_from_file(self, fn, debug_break=False, top_level=False, only_first_frame=False):
        """Create a trajectory based on give filename.

        Args:

           top_level:
                top_level = True, don't make a huge tree of objects (Residues/Atoms) == amazing speed up!
                Useful if you need only frames, energies and coords as text.
                You only get the info that is in header of each frame.

                top_level = False, makes huge tree of objects (Residues/Atoms) == very slow for a huge trajectories

        .. warning:: Loads up whole trafl file into memory, and get stuck. Use this if you want to compute e.g. distances between atoms, get the positions of specified atoms etc. If you can not process your trajectory
        use top_level=True or look at load_from_string() to load a frame by frame from a file.

        h(eader), l(line), f(ile)
        """
        self.frames = []
        f = (line for line in open(fn))
        h = next(f).strip()
        l = next(f).strip()
        c = 0
        self.frames.append(Frame(c, h, l, top_level))

        if not only_first_frame:
            while 1:
                c += 1
                try:
                    h = next(f).strip()
                except StopIteration:
                    break
                l = next(f).strip()
                if h and l:
                    self.frames.append(Frame(c, h, l, top_level))
                    if debug_break:
                        break
                    if c % 1000 == 0:
                        print(c / 1000, 'k loaded...')
                        gc.collect()

    def load_from_string(self, c, txt):
        """Create a trajectory based on given string (txt) with id given by c.

        Faster method, loads only one frame at a time to memory, and after computations
        loads the next frame (memory efficient)."""
        h, l = txt.split('\n')
        self.frames.append(Frame(c, h, l))

    def sort(self, inplace=True):
        """Sort frames within the trajectory according to energy."""
        def getEnergy(frame):
            return frame.energy
        frames_sorted = sorted(self.frames, key=getEnergy)
        if inplace:
            self.frames = frames_sorted
        else:
            return frames_sorted

    def load_from_list(self, frames):
        self.frames = frames

    def save(self, fn, verbose=True):
        """Save the trajectory to file."""
        with open(fn, 'w') as fi:
            for f in self.frames:
                # [295.0, 5.0, -1428.683789, -1435.583554, 0.9]
                fi.write(' '.join([str(x) for x in f.header]) + '\n')
                fi.write(f.coords + '\n')
        if verbose:
            print('Saved to ' + fn)

    def plot_energy(self, plotfn='plot.png'):
        """Plots the SimRNA energy of the trajectory.

        .. image:: ../pngs/simrnatrajectory.png
        """

        # plotting inside ipython
        import matplotlib.pyplot as plt
        import matplotlib
        plt.plot([f.energy for f in self.frames])
        plt.ylabel('# frames')
        plt.ylabel('energies')
        plt.title('SimRNATrajectory: energies over frames')
        plt.grid()
        plt.savefig(plotfn, figsize=(30, 10))  # bbox_inches='tight',
        plt.close()

    def __len__(self):
        """Get number of frames"""
        return len(self.frames)

    def __getitem__(self, i):
        return self.frames[i]


class Frame:
    """Frame

    Syntax of header:
     - write_number
     - replica_id
     - total_energy
     - energy_without_restraints
     - temperature

    .. warning:: If there is an invalid frame, please use `repair_trafl.py` to fix the trajectory first.
    """

    def __init__(self, id, header, coords, top_level=False):
        """
        top_level = False, makes huge tree of objects (Residues/Atoms) == slower
        Important when you want to compute e.g. disdances between atoms, get the positions of specified atoms etc.
        """
        self.id = id
        self.header = [float(h) for h in header.split()]
        # 443 10 -1211.659921 -1218.294753 1.100000
        l = header.split()
        if len(l) != 5:
            raise Exception('Invalid frame, please use `repair_trafl.py` to fix it.')
        self.energy = float(l[3])

        self.coords = coords
        coords = deque(coords.split())
        if top_level:
            self.residues = []
        else:
            self.residues = []
            c = 1
            while coords:
                p = Atom('p', coords.popleft(), coords.popleft(), coords.popleft())
                c4p = Atom('c4p', coords.popleft(), coords.popleft(), coords.popleft())
                n1n9 = Atom('n1n9', coords.popleft(), coords.popleft(), coords.popleft())
                b1 = Atom('b1', coords.popleft(), coords.popleft(), coords.popleft())
                b2 = Atom('b2', coords.popleft(), coords.popleft(), coords.popleft())
                r = Residue(c, p, c4p, n1n9, b1, b2)
                self.residues.append(r)
                c += 1

    def __repr__(self):
        return 'Frame #' + str(self.id) + ' e:' + str(round(self.energy, 2))

    def __len__(self):
        """Get a number of residues"""
        return len(self.residues)


class Residue:
    """Create Residue object.

       Each residue in SimRNA coarse-grained represantation consists only 5 coarse-grained atoms:

       - backbone: p = phospate group, c4p = sugar moiety
       - nucleotide: n1n9 = N1 for pyrimidines, N9 for purines, b1 = C2 for purines and pyrimidines, \
       b2 = C4 for pyrimidines, C6 for purines
    """

    def __init__(self, id, p, c4p, n1n9, b1, b2):
        self.id = id
        self.p = p
        self.c4p = c4p
        self.n1n9 = n1n9
        self.b1 = b1
        self.b2 = b2
        self.atoms = [p, c4p, n1n9, b1, b2]
        # self.center

    def __sub__(self, other_residue):
        diff = self.get_center() - other_residue.get_center()
        return np.sqrt(np.dot(diff, diff))

    def get_atoms(self):
        """Return all atoms"""
        return self.atoms

    def get_center(self):
        """Return MB for residue ```((self.n1n9 + self.b2) / 2)```"""
        return (self.n1n9 + self.b2) / 2

    def __repr__(self):
        return 'Residue ' + str(self.id)


class Atom:
    """Atom
    x
    y
    z
    coord
    """

    def __init__(self, name, x, y, z):
        """Create Atom object."""
        self.name = name
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.coord = np.array([self.x, self.y, self.z])

    def get_coord(self):
        """Return coords (np.array)."""
        return self.coord

    def __repr__(self):
        return 'Atom ' + str(self.x) + ' ' + str(self.y) + ' ' + str(self.z)

    def __sub__(self, other_atom):
        """Calculate distance between two atoms.

        Example::

            distance = atom1 - atom2

        """
        diff = self.coord - other_atom.coord
        return np.sqrt(np.dot(diff, diff))

    def __add__(self, other_atom):
        return self.coord + other_atom.coord


# main
if __name__ == '__main__':
    s = SimRNATrajectory()
    s.load_from_file('test_data/6c2ca958-d3aa-43b2-9f02-c42d07a6f7e9_ALL.trafl', top_level=True)
    s.plot_energy('plot.png')

    s = SimRNATrajectory()
    s.load_from_file(
        'test_data/8b2c1278-ee2f-4ca2-bf4a-114ec7151afc_ALL_thrs6.20A_clust01.trafl', top_level=False)
    for f in s.frames:
        print('header of each frame:', f.header)
        # prints out the energy of restraints (total_energy - energy_without_restraints)
        print(f.header[3] - f.header[4])
        print(f.coords)  # prints coords of each coarse-grained atom

    s2 = SimRNATrajectory()
    traj = """1 1 1252.257530 1252.257530 0.950000
 53.570 23.268 39.971 55.119 24.697 43.283 55.145 27.966 42.270 55.258 29.321 42.618 54.313 29.909 40.572 57.246 41.229 41.492 57.056 39.572 45.104 55.672 36.799 43.722 55.491 33.284 44.069 55.013 33.922 41.769"""

    s2.load_from_string(0, traj)
    for f in s2.frames:
        print(f)

        # print f.coords
        r = f.residues[0]
        print(r.get_atoms())
        print(r.p)  # prints the position of P group

        print(r.c4p)
        print(r.c4p.get_coord())
        print(r.b1 - r.p)  # prints the distance between C4 and P coarse grained-atoms

        print(r.n1n9)
        print(r.b2)
        print(r.get_center())
