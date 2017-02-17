#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SimRNATrajectory module

SimRNATrajectory / Frame / Residue / Atom"""

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
        to load the data."""
        self.frames = []
        
    def load_from_file(self, fn, debug_break=False, top_level=False):
        """Create a trajectory based on give filename.

        top_level = True, don't make a huge tree of objects (Residues/Atoms) == amazing speed up! 
        Useful if you need only frames, energies and coords as text.

        h(eader), l(line), f(ile).
        """
        self.frames = []
        f = (line for line in open(fn).xreadlines())
        h = f.next().strip()
        l = f.next().strip()
        c = 0
        self.frames.append(Frame(c, h, l, top_level))
        while 1:
            c += 1
            try:
                h = f.next().strip()
            except StopIteration:
                break
            l = f.next().strip()
            if h and l:
                self.frames.append(Frame(c, h, l, top_level))
                if debug_break: break
                if c % 1000 == 0:
                    print c/1000,'k loaded...'
                    gc.collect()

    def load_from_string(self, c, txt):
        """Create a trajectory based on given string (txt) with id given by c"""
        h,l = txt.split('\n')
        self.frames.append(Frame(c, h, l))        

    def sort(self):
        """Sort frames within the trajectory according to energy."""
        def getEnergy(frame):
            return frame.energy
        frames_sorted = sorted(self.frames, key=getEnergy)
        return frames_sorted

    def load_from_list(self, frames):
        self.frames = frames
        
    def save(self, fn, verbose=True):
        """Save the trajectory to file."""
        with open(fn, 'w') as fi:
            for f in self.frames:
                fi.write(f.header + '\n')
                fi.write(f.coords + '\n')
        if verbose: print('Saved to ' + fn)
            
    def plot_energy(self, plotfn):
        """
        .. image:: ../pngs/simrnatrajectory.png
        """
        
        #plotting inside ipython
        import matplotlib.pyplot as plt
        import matplotlib

        plt.plot([f.energy for f in self.frames])
        plt.ylabel('# frames')
        plt.ylabel('energies')
        plt.title('SimRNATrajectory: energies over frames')
        plt.grid()
        plt.savefig(plotfn,figsize=(30,10)) #  bbox_inches='tight', 

class Frame:
    """Frame

    Syntax of header: 
     - write_number
     - replica_id
     - total_energy 
     - energy_without_restraints
     - temperature
    """
    def __init__(self, id, header, coords, top_level=True):
        self.id = id
        self.header = header
        self.energy = float(header.split(' ')[3])
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
        return 'Frame #' + str(self.id) + ' e:' + str(round(self.energy,2))
                
class Residue:
    """Create Residue object.
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
        return 'Atom ' +  str(self.x) + ' ' + str(self.y) + ' '+ str(self.z)
        
    def __sub__(self, other_atom):
        """Calculate distance between two atoms.

        Example::

            >>> distance=atom1-atom2

        """

        diff = self.coord - other_atom.coord
        return np.sqrt(np.dot(diff, diff)) 

    def __add__(self, other_atom):
        return self.coord + other_atom.coord

#main
if __name__ == '__main__':
    s = SimRNATrajectory()
    s.load_from_file('test_data/6c2ca958-d3aa-43b2-9f02-c42d07a6f7e9_ALL.trafl', top_level=True)
    s.plot_energy('plot.png')

    s = SimRNATrajectory()
    s.load_from_file('test_data/8b2c1278-ee2f-4ca2-bf4a-114ec7151afc_ALL_thrs6.20A_clust01.trafl')
    for f in s.frames:
        #print f.header
        #print f.coords
        r = f.residues[0]
        print r.get_atoms()
        print r.p
        print r.c4p
        print r.c4p.get_coord()
        print r.c4p - r.p
        print r.n1n9
        print r.b2
        print r.get_center()
        break

    s2 = SimRNATrajectory()
    traj = """1 1 1252.257530 1252.257530 0.950000
 53.570 23.268 39.971 55.119 24.697 43.283 55.145 27.966 42.270 55.258 29.321 42.618 54.313 29.909 40.572 57.246 41.229 41.492 57.056 39.572 45.104 55.672 36.799 43.722 55.491 33.284 44.069 55.013 33.922 41.769"""
    s2.load_from_string(0, traj)
    for f in s2.frames:
        print f
