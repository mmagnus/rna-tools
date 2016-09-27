#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SimRNATrajectory module"""

from collections import deque
import numpy as np

import gc


class SimRNATrajectory:
    """SimRNATrajectory:

    -> Frame / Residue / Atom"""
    def __init__(self):
        self.frames = []
        
    def load_from_file(self, fn):
        """h(eader), l(line), f(ile)"""
        frames = []
        f = (line for line in open(fn).xreadlines())
        c = 1
        h = f.next().strip()
        l = f.next().strip()
        frames.append(Frame(c, h, l))
        c = 0
        while 1:
            c += 1
            try:
                h = f.next().strip()
            except StopIteration:
                break
            l = f.next().strip()
            if h and l:
                frames.append(Frame(c, h, l))
                if c % 1000 == 0:
                    print c/1000,'k cleaning...'
                    gc.collect()
        self.frames = frames
        #with open(fn) as f:
        #    for line in f:
        #        do_something_with(line)
        #f = open(fn)  # how much is loaded in memory?
        #f.next()  

    def load_from_string(self, txt):
        c = 0
        h,l = txt.split('\n')
        self.frames.append(Frame(c, h, l))        
    
class Atom:
    """Atom"""
    def __init__(self, name, x, y, z):
        self.name = name
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.coord = np.array([self.x, self.y, self.z])

    def __repr__(self):
        return 'Atom ' +  str(self.x) + ' ' + str(self.y) + ' '+ str(self.z)

    def get_coord(self):
        return self.coord
        
    def __sub__(self, other_atom):
        diff = self.coord - other_atom.coord
        return np.sqrt(np.dot(diff, diff)) 

    def __add__(self, other_atom):
        return self.coord + other_atom.coord

class Residue:
    """Residue"""
    def __init__(self, id, p, c4p, n1n9, b1, b2):
        self.id = id
        self.p = p
        self.c4p = c4p
        self.n1n9 = n1n9
        self.b1 = b1
        self.b2 = b2
        self.atoms = [p, c4p, n1n9, b1, b2]
        # self.center

    def __repr__(self):
        return 'Residue ' + str(self.id)

    def get_atoms(self):
        return self.atoms

    def get_center(self):
        return (self.n1n9 + self.b2) / 2

class Frame:
    """Frame"""
    def __init__(self, id, header, coords):
        """header: write_number replica_id total_energy energy_without_restraints temperature"""
        self.id = id
        self.header = header
        self.energy = float(header.split(' ')[3])
        self.coords = coords
        coords = deque(coords.split())
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
                
#main
if __name__ == '__main__':
    s = SimRNATrajectory()
    s.load_from_file('8b2c1278-ee2f-4ca2-bf4a-114ec7151afc_ALL_thrs6.20A_clust01.trafl')
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
    s2.load_from_string(traj)
    for f in s2.frames:
        print f
