#!/usr/bin/env python
"""
Cite & authors:

[1]	T. Waleń, G. Chojnowski, P. Gierski, and J. M. Bujnicki, “ClaRNA: a classifier of contacts in RNA 3D structures based on a comparative analysis of various classification schemes.,” Nucleic Acids Research, vol. 42, no. 19, pp. e151–e151, Oct. 2014.
"""
import sys
import os
import re
import math
import gzip
import simplejson as json
import networkx as nx
import itertools
from networkx.readwrite import json_graph
from optparse import OptionParser
from multiprocessing import Pool, Queue, Process

from scipy.spatial import cKDTree

import numpy as np
from numpy import array, dot, sqrt
from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer

from rna_tools.tools.clarna_play.ClaRNAlib.distances import range_value_distance, rmsd_distance, doublet_params_dict, center_vector, bp_distance, vector_length, Residue, Doublet
from rna_tools.tools.clarna_play.ClaRNAlib.structure_ciach import StructureCiachCiach
from rna_tools.tools.clarna_play.ClaRNAlib.utils import *
from rna_tools.tools.clarna_play.ClaRNAlib.cl_settings_global import *

MAX_RES_DIST = 7.0
MAX_C1P_DIST = 14
N_TYPES = ('A','C','G','U')

descriptions_dict = {}

class KeyboardInterruptError(Exception):
    """dummy exception for Ctrl-C"""
    pass

class ContactProgram(object):

    def load_from_json(self, json_fn):
        prg = load_json(json_fn)
        self.commands = self._load_program(prg)

    def _load_expression(self, p):
        op_code = p[0]
        if op_code == "IF":
            assert(len(p)==3)
            cond_prg = self._load_program(p[1])
            then_prg = self._load_program(p[2])
            return MyExprIf(cond_prg, then_prg)
        elif op_code == "ALIGN_TO":
            assert(len(p)==3 or len(p)==4)
            if len(p)==4:
                return MyExprAlign(p[1], p[2], p[2], p[3])
            else:
                return MyExprAlign(p[1], p[2], p[2])
        elif op_code == "ALIGN_TO2":
            assert(len(p)==4 or len(p)==5)
            if len(p)==5:
                return MyExprAlign(p[1], p[2], p[3], p[4])
            else:
                return MyExprAlign(p[1], p[2], p[3])
        elif op_code=="ALIGN_TO_WITH_RMSD_ATOMS" or op_code=="ALIGN_TO_WITH_MULTIPLE_RMSD_ATOMS":
            assert(len(p)==5 or len(p)==6)
            multiple = False
            if op_code=="ALIGN_TO_WITH_MULTIPLE_RMSD_ATOMS":
                multiple = True
            if len(p)==6:
                return MyExprAlign(p[1], p[2], p[3], p[4], rmsd_atoms=p[5], multiple=multiple)
            else:
                return MyExprAlign(p[1], p[2], p[3], rmsd_atoms=p[4], multiple=multiple)
        elif op_code == "COMPUTE_BP_DISTANCE":
            assert(len(p)==2 or len(p)==3)
            return MyExprComputeBpDistance(*p[1:])
        elif op_code == "COMPARE_DISTANCES":
            assert(len(p)==2 or len(p)==3)
            if len(p)==3:
                return MyExprCompareDistances(p[1], p[2])
            else:
                return MyExprCompareDistances(p[1])
        elif op_code == "DISTANCE_TO_POINTS":
            assert(len(p)==3)
            return MyExprDistanceToPoints(p[1],p[2])
        elif op_code == "COMPUTE_DISTANCE_TO_POINTS":
            assert(len(p)==3)
            return MyExprComputeDistanceToPoints(p[1],p[2])
        elif op_code == "LE":
            assert(len(p)==3)
            if isinstance(p[1],list):
                p[1] = self._load_expression(p[1])
            if isinstance(p[2],list):
                p[2] = self._load_expression(p[2])
            return MyExprLE(p[1], p[2])
        elif op_code == "LT":
            assert(len(p)==3)
            return MyExprLT(p[1], p[2])
        elif op_code == "GE":
            assert(len(p)==3)
            return MyExprLE(p[2], p[1])
        elif op_code == "GT":
            assert(len(p)==3)
            return MyExprLT(p[2], p[1])
        elif op_code == "EQ":
            assert(len(p)==3)
            return MyExprEQ(p[1], p[2])
        elif op_code == "OR":
            assert(len(p)>=2)
            return MyExprOr(*[self._load_expression(x) for x in p[1:]])
        elif op_code == "IN_RANGE":
            assert(len(p)==4)
            return MyExprInRange(p[1], p[2], p[3])
        elif op_code == "RETURN":
            assert(len(p)==2 or len(p)==3)
            if len(p)==2:
                return MyExprReturn(p[1])
            elif len(p)==3:
                return MyExprReturn(p[1],p[2])
            else:
                raise Exception("Bad return")
        elif op_code == "RETURN_IF":
            assert(len(p)==5)
            if len(p)==5:
                return MyExprReturnIf(p[1],p[2],p[3],p[4])
            else:
                raise Exception("Bad return if")
        else:
            raise Exception("Unknown expression: %s" % op_code)

    def _load_program(self,prg):
        res = []
        for p in prg:
            res.append(self._load_expression(p))
        return res

    def value(self, state):
        """TODO -- duplikacja kodu!"""
        res = None
        new_state = state
        for c in self.commands:
            (r,new_state) = c.value(new_state)
            if r is not None:
                if res is None:
                    res = []
                if isinstance(r,list):
                    res += r
                else:
                    res += [r]
        return (res, new_state)


class ContactProgramState(object):

    def __init__(self, r1, r2, n_type, options=None, doublet=None):
        self.rmsd = None
        self.bp_distance = None
        self.distance_to_points = None
        self.dist_score = None
        self.reference_doublet_id = None
        self.doublet = None

        if doublet:
            self.doublet = doublet
            self.n_type = doublet.n_type
            # dodac do parametrow!
            self.consecutive = self.is_consecutive(doublet.res1.points, doublet.res2.points)
            self.simple_n_type = self._simplify_n_type(doublet.n_type)
            ######################
            self.r1 = doublet.res1.points
            self.r2 = doublet.res2.points
        else:
            self.r1 = r1
            self.r2 = r2
            self.n_type = n_type
            self.simple_n_type = self._simplify_n_type(n_type)
            self.consecutive = self.is_consecutive(r1,r2)

            if re.match('^[ACGU]{2}$',n_type):
                self.normalized_r1, self.normalized_r2 = normalize_points((self.r1,self.r2), n_type[0])
                if self.normalized_r1 is not None:
                    self.normalized_r1['center'] = center_vector(self.normalized_r1, n_type[0])
                if self.normalized_r2 is not None:
                    self.normalized_r2['center'] = center_vector(self.normalized_r2, n_type[1])
            else:
                self.normalized_r1 = {}
                self.normalized_r2 = {}
            self.doublet_params = doublet_params_dict((self.r1,self.r2), n_type)
            self.doublet_params_ph = doublet_params_dict((self.r1,self.r2), n_type, 'base-phosphate')
            self.doublet_params_br = doublet_params_dict((self.r1,self.r2), n_type, 'base-ribose')

        self.options = options
        if self.options.save_scores or self.options.show_scores_for:
            self.save_all_scores = True
        else:
            self.save_all_scores = False

    @property
    def dist_matrix(self):
        if not hasattr(self,'_dist_matrix'):
             self._dist_matrix = self._compute_distances(list(self.r1.keys()), list(self.r2.keys()))
        return self._dist_matrix

    def is_consecutive(self,r1,r2):
        for a1,a2 in [(r1.get('P'),r2.get("O3'")), (r2.get('P'),r1.get("O3'"))]:
            if a1 is not None and a2 is not None:
                if vector_length(np.array(a1,'f')-np.array(a2,'f'))<1.7:
                    return 1
        return 0
                    

    def _simplify_n_type(self,n_type):
        PURINE = ['A','G']
        PYRIMIDINE = ['C','U']
        res = ""
        for c in n_type:
            if res != "":
                res += "-"
            if c in PURINE:
                res += "Pu"
            elif c in PYRIMIDINE:
                res += "Py"
        return res

    def get_attr(self, a):
        res = self._get_attr(a)

        if False:
            tmp = self.doublet
            self.doublet = None
            res2 = self._get_attr(a)
            self.doublet = tmp

        #if res!=res2:
        #    print "a=%s res=%s res2=%s" % (a,res,res2)

        return res
         
    def _get_attr(self, a):
        if a=="n_type":
            return self.n_type
        elif a=="simple_n_type":
            return self.simple_n_type
        elif a=="rmsd":
            return self.rmsd
        elif a=="bp_distance":
            return self.bp_distance
        elif a=="distance_to_points":
            return self.distance_to_points
        elif a=="dist_score":
            return self.dist_score
        elif a in ["consecutive","is_consecutive"]:
            return self.consecutive
        elif re.match('^d_param_br_.*',a):
            key = a[11:]
            if self.doublet:
                if re.match("^(dist_[A-Z]|i_|ii_|oxygens)",key):
                    if self.doublet.br_info is None:
                        return None
                    return self.doublet.br_info.get(key)
                else:
                    return getattr(self.doublet, key)
            else:
                if self.doublet_params_br is None:
                    return None
                return self.doublet_params_br.get(key)
        elif re.match('^d_param_ph_.*',a):
            key = a[11:]
            if self.doublet:
                if re.match("^(dist_[A-Z]|i_|ii_|oxygens)",key):
                    if self.doublet.ph_info is None:
                        return None
                    return self.doublet.ph_info.get(key)
                else:
                    return getattr(self.doublet, key)
            else:
                if self.doublet_params_ph is None:
                    return None
                return self.doublet_params_ph.get(key)
        elif re.match('^d_param_.*',a):
            if self.doublet:
                return getattr(self.doublet, a[8:])
            else:
                if self.doublet_params is None:
                    return None
                return self.doublet_params.get(a[8:])
        elif re.match('^norm1_.*',a):
            if self.doublet:
                if self.doublet.normalized1 is None:
                    return None
                if re.match("^norm1_center",a):
                    return self.doublet.normalized1.center
                else:
                    return self.doublet.normalized1.points.get(a[6:])
            else:
                if self.normalized_r1 is None:
                    return None
                return self.normalized_r1.get(a[6:])
        elif re.match('^norm2_.*',a):
            if self.doublet:
                if self.doublet.normalized2 is None:
                    return None
                if re.match("^norm2_center",a):
                    return self.doublet.normalized2.center
                else:
                    return self.doublet.normalized2.points.get(a[6:])
            else:
                if self.normalized_r2 is None:
                    return None
                return self.normalized_r2.get(a[6:])
        return None

    def align(self, points, use_atoms_chain_a, use_atoms_chain_b, reference_doublet_id=None, rmsd_atoms=None, multiple=False):
        self.rmsd = 1000.0
        self.reference_doublet_id = reference_doublet_id

        if self.options and self.options.disable_align_to:
            return self
        
        self.rmsd = rmsd_distance((self.r1, self.r2), points, (use_atoms_chain_a, use_atoms_chain_b), rmsd_atoms, multiple)

        return self

    def calc_bp_distance(self, points, reference_doublet_id=None):
        self.bp_distance = 1000.0
        self.reference_doublet_id = reference_doublet_id
        if self.doublet:
            self.bp_distance = bp_distance((self.doublet.normalized1.points, self.doublet.normalized2.points), points, self.n_type, already_normalized=True)
        else:
            self.bp_distance = bp_distance((self.normalized_r1, self.normalized_r2), points, self.n_type, already_normalized=True)
        return self

    def _compute_distances(self,keys1,keys2):
        res = {}
        for k1 in keys1:
            row = {}
            for k2 in keys2:
                row[k2] = dist(self.r1[k1], self.r2[k2])
            res[k1] = row
        return res

    def compare_dist_matrix(self, dist_matrix, reference_doublet_id):
        self.dist_score = 1000.0
        self.reference_doublet_id = reference_doublet_id

        if self.options and self.options.disable_compare_distances:
            return self

        keys1 = list(set(dist_matrix.keys()).intersection(list(self.r1.keys())))
        if len(keys1)==0:
            return self
        keys2 = list(set(dist_matrix[keys1[0]].keys()).intersection(list(self.r2.keys())))
        if len(keys2)==0:
            return self

        curr_distances = self.dist_matrix

        score = 0
        n = 0
        for k1 in keys1:
            for k2 in keys2:
                s = range_value_distance(dist_matrix[k1][k2], curr_distances[k1][k2])
                # print "k1=%s k2=%s s=%.5f dm=%s d=%.5f" % (k1,k2,s,dist_matrix[k1][k2], curr_distances[k1][k2])
                score += s
                n += 1
        score /= n
        self.dist_score = (1-score)
        return self

class EmptyContactProgramState(ContactProgramState):

    def __init__(self,**kwargs):
        for k,v in list(kwargs.items()):
            setattr(self,k,v)
        for k,default_value in (('save_all_scores',False),('rmsd',None),('reference_doublet_id',None),('dist_score',None),):
            if not hasattr(self,k):
                setattr(self,k,default_value)

class MyResult(object):

    def __init__(self, desc, score=None, real_score=None, reference_id=None, ret_score_id=None):
        self.desc = desc
        if score is None:
            score=1.0
        if real_score is None:
            real_score=1000.0
        self.score = score
        self.real_score = real_score
        self.reference_id = reference_id
        self.ret_score_id = ret_score_id

    def __str__(self):
        return "%s [score=%.3f/%.3f,ref=%s,ret_score_id=%s]" % (self.desc, self.score, self.real_score, self.reference_id, self.ret_score_id)

    def __repr__(self):
        return self.__str__()


class MyExpression(object):

    def value(self, state):
        return (None, state)

    def arg_value(self, a, state):
        if isinstance(a,str) and len(a)>1:
            if a[0]=='$':
                return state.get_attr(a[1:])
        elif isinstance(a,MyExpression):
            (r,_state) = a.value(state)
            # TODO we ignore state here!
            return r
        return a

class MyExprReturnIf(MyExpression):

    def __init__(self, test_expr, bound_for_ok, bound_for_fuzzy, ret_value):
        self.test_expr = test_expr
        self.bound_for_ok = bound_for_ok
        self.bound_for_fuzzy = bound_for_fuzzy
        self.ret_value = ret_value

    def value(self, state):
        ret_reference = None
        if state.reference_doublet_id is not None:
            ret_reference = state.reference_doublet_id

        expr = MyExpression()
        ret_score = expr.arg_value(self.test_expr, state)
        if ret_score <= self.bound_for_ok:
            r = (1-0.5*ret_score/max(0.00000001,self.bound_for_ok))
            r = max(0,min(1,r))
            return (MyResult(self.ret_value, r, ret_score, ret_reference, self.test_expr), state)
        if ret_score <= self.bound_for_fuzzy:
            bound_diff = max(0.1, self.bound_for_fuzzy - self.bound_for_ok)
            return (MyResult('?'+self.ret_value, 0.5-0.5*(ret_score-self.bound_for_ok)/bound_diff, ret_score, ret_reference, self.test_expr), state)
        if state.save_all_scores:
            return (MyResult('!'+self.ret_value, 0, ret_score, ret_reference, self.test_expr), state)
        return (None, state)

class MyExprReturn(MyExpression):

    def __init__(self, ret_value, ret_score=None):
        self.ret_value = ret_value
        self.ret_score = ret_score

    def value(self, state):
        ret_desc = self.ret_value
        ret_score = None
        if self.ret_score is not None:
            expr = MyExpression()
            ret_score =  expr.arg_value(self.ret_score, state)
        ret_reference = None
        if state.reference_doublet_id is not None:
            ret_reference = state.reference_doublet_id

        result = MyResult(ret_desc, ret_score, ret_score, ret_reference, self.ret_score)
        return (result,state)

class MyExprAlign(MyExpression):

    def __init__(self, points, use_atoms_chain_a, use_atoms_chain_b, reference_doublet_id=None, rmsd_atoms=None, multiple=False):
        self.points = points
        self.use_atoms_chain_a = use_atoms_chain_a
        self.use_atoms_chain_b = use_atoms_chain_b
        self.reference_doublet_id = reference_doublet_id
        self.rmsd_atoms = rmsd_atoms
        self.multiple = multiple

    def value(self, state):
        new_state = state.align(self.points, self.use_atoms_chain_a, self.use_atoms_chain_b, self.reference_doublet_id, self.rmsd_atoms, self.multiple)
        return (None, new_state)

class MyExprComputeBpDistance(MyExpression):

    def __init__(self, points, reference_doublet_id=None):
        self.points = points
        self.reference_doublet_id = reference_doublet_id

    def value(self, state):
        new_state = state.calc_bp_distance(self.points, self.reference_doublet_id)
        return (None, new_state)

class MyExprCompareDistances(MyExpression):

    def __init__(self, dist_matrix, reference_doublet_id=None):
        self.dist_matrix = dist_matrix
        self.reference_doublet_id = reference_doublet_id

    def value(self, state):
        new_state = state.compare_dist_matrix(self.dist_matrix, self.reference_doublet_id)
        return (None, new_state)


class MyExprIf(MyExpression):

    def __init__(self, conditions, then_prg):
        self.conditions = conditions
        self.then_prg = then_prg

    def value(self, state):
        cond_value = True
        new_state = state
        for c in self.conditions:
            (c_res, new_state) = c.value(new_state)
            if not c_res:
                cond_value = False
                break
        if cond_value:
            res = None
            for c in self.then_prg:
                (r,new_state) = c.value(new_state)
                if r is not None:
                    if res is None:
                        res = []
                    if isinstance(r,list):
                        res += r
                    else:
                        res += [r]
            return (res, new_state)
        else:
            return (None, new_state)

class MyExprOr(MyExpression):

    def __init__(self, *conditions):
        self.conditions = conditions

    def value(self, state):
        res = False
        new_state = state
        for c in self.conditions:
            (c_res, new_state) = c.value(new_state)
            if c_res:
                res = True
                break
        return (res, new_state)

class MyExpr2arg(MyExpression):

    def __init__(self, arg1, arg2):
        self.arg1 = arg1
        self.arg2 = arg2

class MyExpr3arg(MyExpression):

    def __init__(self, arg1, arg2, arg3):
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3


class MyExprLE(MyExpr2arg):

    def value(self, state):
        v1 = self.arg_value(self.arg1,state)
        v2 = self.arg_value(self.arg2,state)
        result = (v1 <= v2)
        return (result, state)

class MyExprDistanceToPoints(MyExpr2arg):

    def __init__(self, arg1, arg2):
        super(MyExprDistanceToPoints, self).__init__(arg1,arg2)
        self.tree = cKDTree(arg2)
        
    def value(self, state):
        result = 999.0
        p = self.arg_value(self.arg1,state)
        if p is not None:
            d,j = self.tree.query(p,k=1)
            result = d
        return (result, state)

class MyExprComputeDistanceToPoints(MyExprDistanceToPoints):

    def value(self, state):
        (result, state) =  super(MyExprComputeDistanceToPoints, self).value(state)
        state.distance_to_points = result
        return (None, state)

class MyExprLT(MyExpr2arg):

    def value(self, state):
        result = (self.arg_value(self.arg1,state) < self.arg_value(self.arg2,state))
        return (result, state)


class MyExprEQ(MyExpr2arg):

    def value(self, state):
        v1 = self.arg_value(self.arg1,state)
        v2 = self.arg_value(self.arg2,state)
        result = (v1==v2)
        return (result, state)

class MyExprInRange(MyExpr3arg):

    def value(self, state):
        a1_val = self.arg_value(self.arg1,state)
        a2_val = self.arg_value(self.arg2,state)
        a3_val = self.arg_value(self.arg3,state)
        if a1_val is None or a2_val is None or a3_val is None:
            return (False, state)
        result = (a1_val >= a2_val and a1_val <= a3_val)
        return (result, state)



################################

def valid_pair_new2(curr_points, points, attrs):
    eps = 0.3
    atoms = ["C1'","C2","C4","C6","N1","N3"]
    atoms1 = ["00-%s"%a for a in atoms]+["10-%s"%a for a in atoms]
    # print atoms1, curr_points.keys(), points.keys()
    for a in atoms1:
        if a not in curr_points or a not in points:
            print("missing atom %s" % a)
            return (0,1000,1000,1000)
    p1 = array([points[k] for k in atoms1],'f')
    p2 = array([curr_points[k] for k in atoms1],'f')
    diff = p1-p2
    sup = SVDSuperimposer()
    sup.set(p1, p2)
    sup.run()
    (rot,tran) = sup.get_rotran()
    pp2 = dot(p2, rot)+tran
    good = 0
    bad = 0
    bad_sum = 0
    total_sum = 0
    for i,diff_v in enumerate(diff):
        d = vector_length(d)
        # max_d = max(0.1,attrs[atoms1[i]]['avg_d'])
        max_d = max(0.2,min(0.8,attrs[atoms1[i]]['max_d']))
        total_sum += d
        if d < max_d+eps:
            good += 1
        else:
            bad += 1
            bad_sum += d-max_d
    return (good, bad, bad_sum, total_sum)

def aggregate_results(res):
    ares = {}
    for r in res:
        if re.match(r'^\!',r.desc):
            continue
        d = r.desc.replace('?',"")
        if d not in ares or r.score > ares[d].score:
            ares[d] = r
    ares = [ares[k] for k in sorted(list(ares.keys()), key=lambda x: (-ares[x].score, ares[x].real_score))]
    return ares

def single_result(res):
    if res is None:
        return None
    a = aggregate_results(res)
    if len(a)==0:
        return None
    return a[0]

def get_status(res, other_results, desc_tool):
    desc_fuzzy = ""
    if res is not None and re.match(r'^\?',res.desc):
        desc_fuzzy = res.desc[1:].split(" ")[0]
        
    if res is None or re.match(r'^\?',res.desc):
        desc = ''
    else:
        desc = res.desc.split(" ")[0]

    r_category = desc_tool.get_desc_category(desc)
    (o_group, o_category, o_desc) = desc_tool.interpret_other_results(other_results)
    if o_group == 'valid':
        if desc == o_desc:
            return 'OK (%s)' % o_category
        elif desc_fuzzy == o_desc:
            return 'OK-FUZZY (%s)' % o_category
        elif desc == '':
            return 'UNDETECTED (%s)' % o_category
        elif r_category != o_category:
            return 'WRONG-GROUP (%s)' % o_category
        else:
            return 'WRONG-DESC (%s)' % o_category
    elif o_group == 'fuzzy':
        if desc in o_desc:
            return 'FUZZY (good-desc-%s)' % r_category
        elif desc=='':
            return 'FUZZY (undetected)'
        elif r_category in o_category:
            return 'FUZZY (good-category)'
        else:
            return 'FUZZY (undetected)'
    elif o_group == 'unclassified':
        if desc!='':
            return 'FALSE-POSITIVE (%s)' % r_category
        else:
            return 'OK (%s)' % r_category

def extract_close_doublets_from_file(name,ids=None,use_old_params=False):
    # TODO: this code should be rewritten, code duplications with StructureCiachCiach
    res_num = {}

    pdb_id = os.path.basename(name).upper()[0:4]
    structure = load_pdb(name)
    ciach = StructureCiachCiach(structure,dont_normalize=True)
    close_doublets = []
    s = []
    for model in structure:        
        s_begin = len(s)
        residues = []
        for chain in model:
            for r in chain:
                rr = simplify_residue(r)
                if rr is not None:
                    next_o3p = ciach._locate_backbone(r,"P")
                    if next_o3p is not None:
                        for a in next_o3p:
                            rr['NEXT:%s'%a.id] = a.get_coord().tolist() 
                    _id = chain.id+str(r.get_id()[1])
                    _num = r.get_id()[1]
                    _resname = r.get_resname().strip()
                    if ids is not None:
                        res_num[_id] = len(s) 
                    if not _resname in ['A','C','G','U']:
                        print("ignoring %s (bad resname: %s)" % (_id,_resname))
                        residues.append(None)
                    else:
                        residues.append(r)
                    if not use_old_params:
                        s.append([chain.id, _resname, _num, Residue(_id,_resname,rr), _id])
                    else:
                        s.append([chain.id, _resname, _num, rr, _id])
                else:
                    print("ignoring %s (missing atoms)" % _id)
        s_end = len(s)
        
        if ids is None:
            for (i,j) in compute_close_doublets(residues, max_distance=MAX_RES_DIST):
                close_doublets.append((i,j,0))
                close_doublets.append((j,i,0))
                
    if ids is not None:
        d_ids = [re.sub("^[^:]+:","",tmp) for tmp in ids]
        for id in d_ids:
            (id1,id2) = id.split(":")
            i = res_num.get(id1)
            j = res_num.get(id2)
            if i is None or j is None:
                print("UNKNOWN doublet %s" % id)
                continue
            close_doublets.append((i,j,0))
    return (pdb_id, s, close_doublets)

#########################

def divide_results_by_cat(res):
    res_by_cat = {'bp':[],'stacking':[],'base-phosphate':[],'base-ribose':[]}
    if res is None:
        return res_by_cat
    for r in res:
        c = DoubletDescTool.get_desc_category(r.desc)
        if c in res_by_cat:
            res_by_cat[c].append(r)
        else:
            if 'unk' not in res_by_cat:
                res_by_cat['unk']=[]
            res_by_cat['unk'].append(r)
    return res_by_cat

#########################

def update_scores_old(scores, n_type, res, other, desc_tool):
    if res is None:
        return
        
    (o_group, o_category, o_desc) = desc_tool.interpret_other_results(other)
    row = {}
    
    for r in res:
        short_desc = re.sub('^[\?\!]','',r.desc.split(" ")[0])+"/"+r.ret_score_id
        if short_desc not in row or row[short_desc][0]>r.real_score:
            row[short_desc]=(r.real_score,r.ret_score_id.replace("$",""))
            
    for short_desc_tmp,(score,score_id) in list(row.items()):
        short_desc = short_desc_tmp.split("/")[0]
        key = short_desc + "/" + n_type + "/" + score_id
        if key not in scores:
            scores[key] = {'valid':[],'possible':[],'invalid':[],'unclassified':[]}
        if o_group=='valid':
            if o_desc==short_desc:
                scores[key]['valid'].append(score)
            else:
                scores[key]['invalid'].append(score)
        elif o_group=='fuzzy':
            if short_desc in o_desc:
                scores[key]['possible'].append(score)
            else:
                scores[key]['invalid'].append(score)
        elif o_group=='unclassified':
            scores[key]['unclassified'].append(score)


def update_scores(scores, doublet_id, n_type, res, groups_tool):
    if res is None:
        return
            
    res_dict = {}
    for r in res:
        short_desc = re.sub('^[\?\!]','',r.desc).split(" ")[0]
        ret_score_id = r.ret_score_id
        if ret_score_id is None:
            ret_score_id = ""
        ret_score_id = ret_score_id.replace("$","")
        key = (short_desc,ret_score_id)
        if key not in res_dict:
            res_dict[key] = []
        res_dict[key].append(r)
    for (short_desc, ret_score_id),res_elems in list(res_dict.items()):
        res_elems.sort(key=lambda x: x.real_score)
        category = DoubletDescTool.get_desc_category(short_desc)
        exp_desc = groups_tool.get_ref_desc(doublet_id, category)
        fuzzy_exp_desc = groups_tool.get_fuzzy_desc(doublet_id, category)
        is_classified = exp_desc!='' or len(fuzzy_exp_desc)>0
        for i,r in enumerate(res_elems):
            key = re.sub('^[\?\!]','',r.desc) + "/" + n_type + "/" + ret_score_id
            if key not in scores:
                scores[key] = {'valid':[],'possible':[],'invalid':[],'unclassified':[]}
            if exp_desc != '' and exp_desc == short_desc:
                if i==0:
                    scores[key]['valid'].append(r.real_score)
            elif short_desc in fuzzy_exp_desc:
                if i==0:
                    scores[key]['possible'].append(r.real_score)
            else:
                if is_classified:
                    scores[key]['invalid'].append(r.real_score)
                else:
                    scores[key]['unclassified'].append(r.real_score)
        

#########################

def normalize_graph(g):

    def better_desc(d1,d2):
        if d1 is not None and d2 is not None:
            if d1[0]!='?' and d2[0]=='?':
                return 1
            elif d1[0]=='?' and d2[0]!='?':
                return 2
            # both not fuzzy
            return 1
        elif d1 is not None and d2 is None:
            return 1
        elif d1 is None and d2 is not None:
            return 2
        else:
            return None
    
    unique_edges = set()
    all_edges = {}
    for (source,target,data) in g.edges(data=True):
        if data['type']!='contact' or data['is_alt']:
            continue
        k = (source,target)
        key = (min(source,target), max(source,target))
        unique_edges.add(key)
        if k not in all_edges:
            all_edges[k] = []
        all_edges[k].append(data)
    for source,target in unique_edges:
        k1 = (source,target)
        k2 = (target,source) 
        edges1 = all_edges.get(k1,[])
        edges2 = all_edges.get(k2,[])
        desc1 = dict([(DoubletDescTool.get_desc_category(e['desc']), e['desc']) for e in edges1])
        edge1 = dict([(DoubletDescTool.get_desc_category(e['desc']), e) for e in edges1])
        desc2 = dict([(DoubletDescTool.get_desc_category(e['desc']), e['desc']) for e in edges2])
        edge2 = dict([(DoubletDescTool.get_desc_category(e['desc']), e) for e in edges2])
        for sc in ['bp','stacking']:
            d1 = desc1.get(sc)
            d2 = desc2.get(sc)
            if d1 is not None:
                rev_d1 = DoubletDescTool.reverse_desc(d1)
            else:
                rev_d1 = None
            if d2 is not None:
                rev_d2 = DoubletDescTool.reverse_desc(d2)
            else:
                rev_d2 = None
            if (d1 is None and d2 is None) or d1==rev_d2:
                continue

            e1 = edge1.get(sc)
            e2 = edge2.get(sc)
                
            b = better_desc(d1,d2)
            print("# conflict %s:%s - d1=%s d2=%s b=%s" % (source,target,d1,d2,b))
            if b==1:
                if d2 is not None:
                    print("NORM: removing contact %s:%s - %s" % (target,source,d2))
                    g.remove_edge(target, source, key=d2)
                print("NORM: adding contact %s:%s - %s" % (target,source,rev_d1))
                g.add_edge(target, source, type='contact', prg='MY', n_type=DoubletDescTool.reverse_n_type(e1['n_type']),
                           desc=rev_d1, full_desc="%s (r)"%rev_d1, key=rev_d1, weight=e1['weight'])
            elif b==2:
                if d1 is not None:
                    print("NORM: removing contact %s:%s - %s" % (source,target,d1))
                    g.remove_edge(source,target, key=d1)
                print("NORM: adding contact %s:%s - %s" % (source,target,rev_d2))
                g.add_edge(source,target, type='contact', prg='MY', n_type=DoubletDescTool.reverse_n_type(e2['n_type']),
                           desc=rev_d2, full_desc="%s (r)"%rev_d2, key=rev_d2, weight=e2['weight'])
    for source in g.nodes():
        out_edges = [(v1,v2,x) for v1,v2,x in g.out_edges([source], data=True) if x.get('type')=='contact' and not re.match('^\?',x.get('desc',''))]
        out_edges = [(v1,v2,x) for v1,v2,x in out_edges if DoubletDescTool.get_desc_category(x.get('desc',''))=='bp']
        if len(out_edges)>1:
            print("NORM: potential conflict on residue %s, got following contacts: %s"  % (source, ["%s:%s=%s"%(v1,v2,x['desc']) for v1,v2,x in out_edges]))
            if False:
                if any([x['desc']=='WW_cis' for v1,v2,x in out_edges]):
                    for (v1,v2,x) in out_edges:
                        if x['desc']!='WW_cis' and re.match('^(W[HS]|SH|SW)_.*',x['desc']):
                            print("NORM: removing contact %s:%s - %s (conflict with WW_cis)" % (v1,v2,x['desc']))
                            g.remove_edge(v1,v2, key=x['desc'])
    return g

def _run_cl(args):
    i,j,l = args
    n_type = "%s%s" % (gl_s[i][1], gl_s[j][1])
    pair_id = "%s%s:%s%s" % (gl_s[i][0],gl_s[i][2], gl_s[j][0],gl_s[j][2])

    if not gl_options.use_old_params:
        (res,_state) = gl_libs[l].value(ContactProgramState(None, None, n_type, gl_options, doublet=Doublet(pair_id,gl_s[i][3], gl_s[j][3])))
    else:
        (res,_state) = gl_libs[l].value(ContactProgramState(gl_s[i][3], gl_s[j][3], n_type, gl_options))
    return res
    
def _find_contacts_using_library(pdb_id, s, close_doublets, libs, options):
    desc_tool = DoubletDescTool(descriptions_dict)
    if options.compare_with:
        groups_tool = GroupsTool(options.compare_with)
    else:
        groups_tool = GroupsTool()

    n = len(s)
    g = nx.MultiDiGraph()
    for i in range(n):
        g.add_node(s[i][4],resname=s[i][1])

    scores = {}
    
    libs_num = len(libs)
    
    if options.threads > 1:
        global gl_s, gl_libs, gl_options
        gl_s = s
        gl_libs = libs
        gl_options = options

        pool = Pool(options.threads)
        try:
            r = [(i,j,l) for (i,j),l in itertools.product([(i,j) for i,j,_d in close_doublets],list(range(libs_num)))]
            _result = pool.map(_run_cl, r)
            result = {}
            for z in range(len(r)):
                result[r[z]] = _result[z]
            _result = None
            pool.close()
        except KeyboardInterrupt:
            print('control-c pressed')
            pool.terminate()
            print('pool is terminated')
        except Exception as e:
            print('got exception: %r, terminating the pool' % (e,))
            traceback.print_exc()
            pool.terminate()
            print('pool is terminated')
        finally:
            pool.join()
            
    for i,j,d in close_doublets:
        n_type = "%s%s" % (s[i][1], s[j][1])
        # n_type = n_type.upper()

        res_i_id = s[i][4]
        res_j_id = s[j][4]
        
        curr_id = pdb_id + ":" +  res_i_id + ":" + res_j_id
        curr_id_rev = pdb_id + ":" + res_j_id + ":" + res_i_id
        
        if options.ignore_bad_doublets:
            if not (s[i][1] in N_TYPES) or not (s[j][1] in N_TYPES):
                print("ignoring %s, invalid n_type %s" % (curr_id, n_type))
                continue
            if re.match('^[zZ]{2}',pdb_id):
                req_atoms = set(['C2','C4','C6'])
            else:
                req_atoms = set(REQ_ATOMS_LIST)
            ok = True
            for rid,ratoms in (res_i_id, s[i][3]),(res_j_id, s[j][3]):
                if not options.use_old_params:
                    missing_atoms = req_atoms.difference(list(ratoms.points.keys()))
                else:
                    missing_atoms = req_atoms.difference(list(ratoms.keys()))
                if len(missing_atoms)!=0:
                    print("ignoring %s, residue %s, missing atoms: %s" % (curr_id, rid, list(missing_atoms)))
                    ok = False
                    break
            if not ok:
                continue                

        g.add_edge(res_i_id, res_j_id, type='dist', weight=d, n_type=n_type)

        pair_id = "%s%s:%s%s" % (s[i][0],s[i][2], s[j][0],s[j][2])
        if options.show_pair and options.show_pair!=pair_id:
            continue

        for l in range(libs_num):
            if options.threads > 1:
                res = result[(i,j,l)]
            else:
                if not options.use_old_params:
                    (res,_state) = libs[l].value(ContactProgramState(None, None, n_type, options, doublet=Doublet(pair_id, s[i][3], s[j][3])))
                else:
                    (res,_state) = libs[l].value(ContactProgramState(s[i][3], s[j][3], n_type, options))
            
            if options.show_scores_for and res is not None:
                new_res = []
                for r in res:
                    if re.match(r'^[\?!]?'+options.show_scores_for, r.desc):
                        new_res.append(r)
                res = new_res
            res_by_cat = divide_results_by_cat(res)
            cur_res_short_dict = {}
            alt_res_short_dict = {}
            for category,results in list(res_by_cat.items()):
                if len(results)>0:
                    cur = single_result(results)
                    cur_res_short_dict[category] = cur
                    def _simple_desc(x):
                        return x.replace("?","").replace("!","").split(" ")[0]
                    d1 = _simple_desc(cur.desc)
                    alt = None
                    for r in sorted(aggregate_results(results), key=lambda x: -x.score):
                        d2 = _simple_desc(r.desc)
                        if d1==d2:
                            continue
                        alt = r
                        break
                    if "?" in cur.desc and alt is not None and alt.score > 0.7*cur.score:
                        alt_res_short_dict[category] = alt

            if options.save_scores:
                if groups_tool.is_valid(curr_id) or groups_tool.is_valid(curr_id_rev):
                    update_scores(scores, curr_id, n_type, res, groups_tool)
                else:
                    print("ignoring %s in updating scores (no info in groups file)" % (curr_id))

            for is_alt,res_dict in (False,cur_res_short_dict),(True,alt_res_short_dict):
                for category,res_short in list(res_dict.items()):
                    if res_short is not None:
                        '''
                        print "%9s, %s: %5s %s-%s %5s %s: %-30s" % (
                                str(i+1)+"_"+str(j+1),
                                s[i][0],
                                s[i][2],
                                s[i][1],
                                s[j][1],
                                s[j][2],
                                s[j][0],
                                res_short
                        )
                        '''
                        d = res_short.desc.split(" ")[0]
                        g.add_edge(res_i_id, res_j_id, type='contact', prg='MY', n_type=n_type,
                                   key=d, desc=d, full_desc=res_short.desc, weight=res_short.score,
                                   reference=res_short.reference_id,is_alt=is_alt)
            if (options.show_pair or options.show_scores_for) and res is not None:
                print("all res: %s" % sorted(res, key=lambda x: (x.score,x.real_score)))
                
    if options.save_scores:
        save_json(options.save_scores, scores, indent=2)
    if options.save_graph:
        if options.normalize_graph:
            g = normalize_graph(g)
        ensure_dir(options.save_graph)
        save_graph(options.save_graph, g, indent=2)
    return g

def find_contacts_using_library_from_doublet_ids(doublet_ids, libs, options):
    ids_by_pdb = {}
    for id in doublet_ids:
        tmp = id.split(":")
        p = tmp[0].upper()
        if p not in ids_by_pdb:
            ids_by_pdb[p]=[]
        ids_by_pdb[p].append(id)
    pdb_id = None
    s = []
    close_doublets = []
    n = 0
    for _pdb_id in sorted(ids_by_pdb.keys()):
        ids = ids_by_pdb[_pdb_id]
        name = os.path.join(options.data_dir, "%s_rna.pdb.gz" % _pdb_id.lower())
        print("loading %s doublets from %s" % (len(ids), name))
        (_p, ss, cc) = extract_close_doublets_from_file(name,ids,use_old_params=options.use_old_params)
        s += ss
        for (i,j,d) in cc:
            close_doublets.append((i+n,j+n,d))
        n = len(s) 
        if pdb_id is None:
            pdb_id = _pdb_id
        else:
            pdb_id = "MIXED"
    return _find_contacts_using_library(pdb_id, s, close_doublets, libs, options)

def find_contacts_using_library_from_file(filename, libs, options):
    print("processing file: %s" % filename)
    bench_start("finding close doublets")
    (pdb_id, s, close_doublets) = extract_close_doublets_from_file(filename,use_old_params=options.use_old_params)
    bench_stop("(found %d close doublets)" % len(close_doublets))
    
    bench_start("find contacts",new_line=True)
    res = _find_contacts_using_library(pdb_id, s, close_doublets, libs, options)
    bench_stop("find contacts")
    return res

def find_contacts_using_library_from_doublet_id(doublet_id, libs, options):
    return find_contacts_using_library_from_doublet_ids([doublet_id], libs, options)

def find_contacts_using_library_from_doublet_ids_fn(doublet_ids_fn, libs, options):
    doublet_ids = load_json(doublet_ids_fn)
    return find_contacts_using_library_from_doublet_ids(doublet_ids, libs, options)

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="""classifier for detecting contacts in RNA
""")
    parser.add_option("--lib", dest="lib",
                  help="read classifier library from file (or comma seprated list of files)", metavar="JSON")
    parser.add_option("-i", "--input", dest="input",
                  help="read input from FILE", metavar="FILE")
    parser.add_option("--input-doublet-id", dest="input_doublet_id",
                  help="process only single doublet", metavar="PDB:R1:R2")
    parser.add_option("--input-doublet-ids", dest="input_doublet_ids",
                  help="process doublets from file", metavar="JSON_FILE")
    parser.add_option("--data-dir", dest="data_dir",
                  help="set the directory with PDB files", metavar="DIR", default='gc-data/pdb_files')
    parser.add_option("--pdb-id", dest="pdb_id",
                  help="set the PDB ids", metavar="PDB_ID")
    parser.add_option("--show-pair", dest="show_pair", metavar="A:B",
                  help="show all results for pair A:B")
    parser.add_option("--normalize-graph", dest="normalize_graph", action='store_true', default=False,
                  help="normalize contact graph")
    parser.add_option("--save-graph", dest="save_graph", metavar="FILE",
                  help="save contact graph")
    parser.add_option("--compare-with", dest="compare_with", metavar="GROUPS_JSON",
                  help="compare with results from other programs")
    parser.add_option("--disable-align-to", dest="disable_align_to", action='store_true',
                  default=False, help="disable align_to commands")
    parser.add_option("--disable-compare-distances", dest="disable_compare_distances", action='store_true',
                  default=False, help="disable compare_distances")
    parser.add_option("--descriptions-dict", dest="descriptions_dict", metavar="JSON_FILE",
                  help="load descriptions dictionary from JSON file")
    parser.add_option("--save-scores", dest="save_scores", metavar="JSON_FILE",
                  help="save scores to JSON file")
    parser.add_option("--show-scores-for", dest="show_scores_for", metavar="DESC",
                  help="show scores for given DESC")
    parser.add_option("--ignore-bad-doublets", dest="ignore_bad_doublets", action='store_true',
                  help="ignore bad doublets")
    parser.add_option("--threads", dest="threads", 
                  help="number of threads", metavar="NUM", default=1)
    parser.add_option("--use-old-params", dest="use_old_params", action='store_true',
                  help="use new params computations", default=False)

    (options, args)  = parser.parse_args()
    if not options.lib:
        lib_dir = os.path.join(os.path.dirname(__file__), "lib")
        cl_groups = ['bp','stacking','base-phosphate','base-ribose','other','other2','other3']
        options.lib = ",".join([os.path.join(lib_dir,"classifier."+x+".json.gz") for x in cl_groups])
    if not options.descriptions_dict:
        options.descriptions_dict = os.path.join(CLARNA_DIR,"descriptions-dict.json")
    options.threads = int(options.threads)
    return (parser, options, args)

def main():
    (parser, options, args) = parse_args()
    if len(args) == 1:
        options.input = args[0]
    if not options.input and not options.input_doublet_id and not options.input_doublet_ids:
        print("specify input")
        parser.print_help()
        exit(1)
    # bench_start("loading descriptions dict")
    global descriptions_dict
    descriptions_dict = load_json(options.descriptions_dict)
    # bench_stop()

    bench_start("loading classifier libraries")
    libs = []
    for lib_fn in options.lib.split(","):
        if not os.path.isfile(lib_fn):
            print("library %s does not exists" % lib_fn)
            parser.print_help()
            exit(1)
        lib = ContactProgram()
        lib.load_from_json(lib_fn)
        libs.append(lib)
    bench_stop()
    if options.input:
        find_contacts_using_library_from_file(options.input, libs, options)
    elif options.input_doublet_id:
        find_contacts_using_library_from_doublet_id(options.input_doublet_id, libs, options)
    elif options.input_doublet_ids:
        find_contacts_using_library_from_doublet_ids_fn(options.input_doublet_ids, libs, options)
    

if __name__ == '__main__':
    main()


