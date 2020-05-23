#!/usr/bin/env python
"""
library for manipulating PDB structures
"""
import sys
import os
import re
import shutil
try:
    from io import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
#import simplejson
import gzip
import tempfile
from optparse import OptionParser
from itertools import combinations
from scipy.spatial import KDTree
from Bio import PDB
import numpy as np
from numpy import array


from rna_tools.tools.clarna_play.ClaRNAlib.normalize_pdb import normalize_pdb
from rna_tools.tools.clarna_play.ClaRNAlib.utils import dist, REQ_ATOMS_LIST, simplify_residue

class StructureCiachCiach(object):

    def __init__(self, structure, dont_normalize=False, req_atoms_list=REQ_ATOMS_LIST):
        self.s = structure
        self.dont_normalize = dont_normalize
        self.res_dict = {}
        # not used at this moment
        # self.backbone_atoms = ["O3'", "C3'", "C4'", "C5'", "O5'", "P"]
        for r in structure.get_residues():
            if r.get_id()[0] != ' ':
                continue
            if dont_normalize:
                key = r.parent.id + (str(r.get_id()[1])+str(r.get_id()[2])).strip()
            else:
                key = str(r.get_id()[1])
            self.res_dict[key] = r
        self.good_residues = []
        self.good_residues_num = {}
        p_coord = []
        o3p_coord = []
        c1p_coord = []
        for (key,r) in list(self.res_dict.items()):
            resname = r.get_resname().strip()
            if len(resname)!=1:
                continue
                
            atoms = simplify_residue(r)
            if atoms is None or not all([a in atoms for a in req_atoms_list]):
                continue

            p_atom = self._locate_atom(r, 'P')
            if p_atom is None:
                p_atom = [0,0,0]
            else:
                p_atom = p_atom.get_coord()
            o3p_atom = self._locate_atom(r, "O3'")
            if o3p_atom is None:
                o3p_atom = [0,0,0]
            else:
                o3p_atom = o3p_atom.get_coord()
            c1p_atom = self._locate_atom(r, "C1'")
            if c1p_atom is None:
                c1p_atom = [0,0,0]
            else:
                c1p_atom = c1p_atom.get_coord()
            
            self.good_residues_num[key] = len(self.good_residues)
            self.good_residues.append(key)
            p_coord.append(p_atom)
            o3p_coord.append(o3p_atom)
            c1p_coord.append(c1p_atom)
        if len(p_coord)==0:
            p_coord.append(array([0,0,0],'f'))
        if len(o3p_coord)==0:
            o3p_coord.append(array([0,0,0],'f'))
        if len(c1p_coord)==0:
            c1p_coord.append(array([0,0,0],'f'))
        self.p_tree = KDTree(p_coord)
        self.o3p_tree = KDTree(o3p_coord)
        self.c1p_coord = c1p_coord
        self.c1p_tree = KDTree(c1p_coord)



    def get_resname(self, num):
        if num not in self.res_dict:
            return '?'
        return self.res_dict[num].get_resname().strip()

    def get_res_atoms_dict(self, num):
        if num not in self.res_dict:
            return None
        res = {}
        for a in self.res_dict[num]:
            res[a.id] = a.get_coord().tolist()
        # add O3' from next nucleotide
        o3p_residue = self._locate_backbone(self.res_dict[num], 'P')
        if o3p_residue is not None:
            for a in o3p_residue:
                res["NEXT:%s"%a.id]=a.get_coord().tolist()
        return res

    def _locate_atom(self, r, id):
        for a in r:
            if a.id==id:
                return a
        return None

    def _locate_backbone(self, r, atom):
        # TODO: use all residues for locating other atoms!
        a = self._locate_atom(r,atom)
        if a is None:
            return None
        point = a.get_coord()
        res_key = None
        if atom=='P':
            (d,p) = self.o3p_tree.query(point, k=1)
        else:
            (d,p) = self.p_tree.query(point, k=1)
        if d<1.4 or d>1.7:
            return None
        res_key = self.good_residues[p]
        rr = self.res_dict[res_key]
        if r.id==rr.id:
            return None
        result = PDB.Residue.Residue(rr.id, rr.resname, rr.segid)
        if atom == 'P':
            other_atoms = ["O3'"]
        else:
            other_atoms = ["P","OP1","O5'","OP2"]

        for a in rr:
            if a.id in other_atoms:
                result.add(a)
        return result

    def initalize_get_neighbours_other():
        self.all_atoms_coords = []
        self.all_atoms_res = []
        for (key,r) in list(self.res_dict.items()):
            for a in r:
                self.all_atoms_coords.append(a.get_coord())
                self.all_atoms_res.append(key)
        self.all_atoms_tree = KDTree(self.all_atoms_coords)

    def residue_distance(self, num1, num2, max_distance=100):
        d = 1000
        if num1 not in self.res_dict or num2 not in self.res_dict:
            return d
        for a1 in self.res_dict[num1]:
           p1 = a1.get_coord()
           for a2 in self.res_dict[num2]:
              p2 = a2.get_coord()
              d = min(d, dist(p1,p2))
              if d<max_distance:
                return d
        return d

    def get_neighbours(self, num, max_distance=4.0):
        if num not in self.good_residues_num:
            return []
        i = self.good_residues_num[num]
        res = []
        points = self.c1p_tree.query_ball_point(self.c1p_coord[i], r=max_distance+15.0)
        for j in points:
            if i==j or j>=len(self.good_residues):
                continue
            num2 = self.good_residues[j]
            if self.residue_distance(num,num2,max_distance)<=max_distance:
                res.append(num2)
        return res

    def get_neighbours_other(self, num, max_distance=4.0):
        if num not in self.res_dict:
            return []
        res = set()
        for a in self.res_dict[num]:
            if not (a.id in ['P',"C1'","C6"]):
                continue
            points = self.all_atoms_tree.query_ball_point(a.get_coord(), r=max_distance)
            for j in set(points):
                if j>=len(self.all_atoms_coords):
                    continue
                num2 = self.all_atoms_res[j]
                if num2!=num:
                    res.add(num2)
        return list(res)


    def get_single_residue(self, num, with_backbone=False):
        if num not in self.res_dict:
            return None

        r1 = self.res_dict[num]

        structure = PDB.Structure.Structure("extracted")
        model = PDB.Model.Model(1)
        structure.add(model)

        c1 = PDB.Chain.Chain("A")
        if with_backbone:
            r = self._locate_backbone(r1,'P')
            if r is not None:
                c1.add(r)
        c1.add(r1)
        if with_backbone:
            r = self._locate_backbone(r1,"O3'")
            if r is not None:
                c1.add(r)
        model.add(c1)
        return structure


    def extract(self,fn, num1, num2, desc, n_type, prg, with_backbone=False):
        if num1 not in self.res_dict:
            return False
        r1 = self.res_dict[num1]
        if num2 not in self.res_dict:
            return False
        r2 = self.res_dict[num2]

        if fn=="dict-no-pdb":
            return {
                'resseq_1': num1[1:],
                'resseq_2': num2[1:],
                'chain_1': num1[0],
                'chain_2': num2[0],
                'desc': desc,
                'n_type': n_type,
                'prg': prg,
            }


        structure = PDB.Structure.Structure("extracted")
        model = PDB.Model.Model(1)
        structure.add(model)

        c1_id = "A"
        c1 = PDB.Chain.Chain(c1_id)
        if with_backbone:
            r = self._locate_backbone(r1,'P')
            if r is not None:
                c1.add(r)
        c1.add(r1)
        if with_backbone:
            r = self._locate_backbone(r1,"O3'")
            if r is not None:
                c1.add(r)

        c2_id = "B"
        c2 = PDB.Chain.Chain(c2_id)
        if with_backbone:
            r = self._locate_backbone(r2,'P')
            if r is not None:
                c2.add(r)
        c2.add(r2)
        if with_backbone:
            r = self._locate_backbone(r2,"O3'")
            if r is not None:
                c2.add(r)

        model.add(c1)
        model.add(c2)

        if fn=="dict":
            output = StringIO.StringIO()
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(output)
            return {
                'resseq_1': num1[1:],
                'resseq_2': num2[1:],
                'chain_1': num1[0],
                'chain_2': num2[0],
                'pdb': output.getvalue(),
                'desc': desc,
                'n_type': n_type,
                'prg': prg,
            }
        else:
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(fn)
            return True
