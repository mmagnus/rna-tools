#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import deque
import numpy as np

class Atom:
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
    def __init__(self, id, header, coords):
        """header: write_number replica_id total_energy energy_without_restraints temperature"""
        print 'header:', header
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
                
class SimRNATrajectory:
    def __init__(self, fn):
        """h(eader), l(line), f(ile)"""
        frames = []
        f = open(fn)
        c = 1
        h = f.readline().strip()
        l = f.readline().strip()
        frames.append(Frame(c, h, l))
        while h:
            c += 1
            h = f.readline().strip()
            l = f.readline().strip()
            if h and l:
                frames.append(Frame(c, h, l))
        print frames
        self.frames = frames
        #with open(fn) as f:
        #    for line in f:
        #        do_something_with(line)
        #f = open(fn)  # how much is loaded in memory?
        #f.next()  
def start(): pass

if __name__ == '__main__':
    s = SimRNATrajectory('8b2c1278-ee2f-4ca2-bf4a-114ec7151afc_ALL_thrs6.20A_clust01.trafl')
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
