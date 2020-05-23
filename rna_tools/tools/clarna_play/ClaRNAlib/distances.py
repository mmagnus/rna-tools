"""library with formulas for calculating the distance scores"""
import math
import re
import itertools
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import scipy

NORMALIZED_BASE = {
  'A': {'N7': [-0.7729614406105143, 2.1201797334590102, 0.0], 'H9': [0.0, -1.0093196646068083, 0.0], '1H6': [1.9209447414524516, 5.150733145313508, 0.0], 'H2': [4.2714470975421825, 1.3179237109436448, 0.0], 'H8': [-2.1004626204274777, 0.44714489628299675, 0.0], '2H6': [0.23043597510468228, 4.699718288132125, 0.0], 'N9': [0.0, 0.0, 0.0], 'C8': [-1.0899103005137063, 0.8335872143873861, 0.0], 'C5': [0.6063291460346252, 2.1057981896577656, 0.0], 'C2': [3.202584142514681, 1.51596183391717, 0.0], 'N1': [2.8648446154295764, 2.824390566297474, 0.0], 'N3': [2.409029202573013, 0.4396032568253272, 0.0], 'C6': [1.563644845822768, 3.1435063393974403, 0.0], 'N6': [1.2031042415497286, 4.447033338301733, 0.0], 'C4': [1.1184904351494842, 0.8031271160408657, 0.0]},
  'C': {'1H4': [-0.5279700289179345, 4.60220182612121, 0.0], 'H1': [0.0, -1.010258706209454, 0.0], 'H6': [-2.082732930782349, 0.12363197089251765, 0.0], '2H4': [1.2038334533450203, 4.411995903837219, 0.0], 'C6': [-1.162053283211799, 0.6983551068687612, 0.0], 'H5': [-2.0367724443934563, 2.647095419696937, 0.0], 'C2': [1.28933438636058, 0.5992529922671946, 0.0], 'N1': [0.0, 0.0, 0.0], 'N3': [1.312136920098135, 1.977774853061991, 0.0], 'N4': [0.28135537444115655, 4.009202043572646, 0.0], 'O2': [2.2694499775195744, -0.12078289494958472, 0.0], 'C5': [-1.1289628161952552, 2.058236947549645, 0.0], 'C4': [0.18336964933177424, 2.654626399935204, 0.0]},
  'U': {'N3': [1.1953166017937347, 1.9738920601305123, 0.0], 'H3': [2.0921523976816854, 2.445178553234257, 0.0], 'H1': [0.0, -1.008947070985391, 0.0], 'H6': [-2.082780056506479, 0.11183564649907562, 0.0], 'C6': [-1.1777603199476516, 0.7096949325723031, 0.0], 'H5': [-2.1245767823028237, 2.616536750689694, 0.0], 'C2': [1.2602866366915804, 0.58925012217277, 0.0], 'N1': [0.0, 0.0, 0.0], 'O4': [0.17715643875791365, 4.020059648793736, 0.0], 'O2': [2.292490396543795, -0.04839432941443833, 0.0], 'C5': [-1.1961299637178628, 2.0625141391615487, 0.0], 'C4': [0.060003385436139034, 2.8083560986778724, 0.0]},
  'G': {'1H2': [4.872516232460849, 0.09778091003631939, 0.0], 'N7': [-0.7513911514240377, 2.1201318166099155, 0.0], 'C2': [3.2594124387538095, 1.333440824418896, 0.0], 'H9': [0.0, -1.0095288732671288, 0.0], 'H1': [3.613335050886427, 3.3745474050771, 0.0], '2H2': [5.295175245866893, 1.7822826455998741, 0.0], 'H8': [-2.101572214696437, 0.4635842689564793, 0.0], 'N9': [0.0, 0.0, 0.0], 'C8': [-1.0886393160814454, 0.841883274438163, 0.0], 'N1': [2.8974071829603005, 2.6580453952999106, 0.0], 'O6': [1.3969862391333745, 4.394213054302783, 0.0], 'N2': [4.595966940484209, 1.0633818952456446, 0.0], 'N3': [2.4038183708450203, 0.34099632542447367, 0.0], 'C6': [1.5671222764357484, 3.191068234453134, 0.0], 'C5': [0.6252146444859401, 2.096901945838508, 0.0], 'C4': [1.117587071308519, 0.7907046767118864, 0.0]},
}

BASE_ATOMS = {
  'A': ['N9', 'C4', 'N3', 'N1', 'C6', 'N6', 'C8', 'C5', 'C2', 'N7'],
  'C': ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C6', 'C5'],
  'U': ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C6', 'C5'],
  'G': ['N9', 'C4', 'N3', 'N1', 'C6', 'O6', 'C8', 'C5', 'C2', 'N7', 'N2'],
}

NORMAL_SUPPORT = {
        'C':['N1','C2','N3','C4','C5','C6'],
        'U':['N1','C2','N3','C4','C5','C6'],
        'G':['N1','C2','C4','N3','C5','C6'],
        'A':['N1','C2','C4','N3','C5','C6'],
        }

PH_OXYGENS = ["OP1","OP2","O5'","NEXT:O3'"]
BR_OXYGENS = ["O2'","O3'","O4'"]

PH_INTERACTIONS = {
    'A': (('C2','H2'), ('C8','H8'), ('N6','1H6'), ('N6','2H6')),
    'C': (('C6','H6'), ('C5','H5'), ('N4','1H4'), ('N4','2H4')),
    'G': (('N1','H1'), ('C8','H8'), ('N2','1H2'), ('N2','2H2')),
    'U': (('C5','H5'), ('N3','H3'), ('C6','H6')),
}
BR_INTERACTIONS = PH_INTERACTIONS

ARCPI = 180.0/np.pi

def normalize_points(points, n_type):
    assert n_type in NORMALIZED_BASE
    assert n_type in BASE_ATOMS
    p1,p2 = points
    norm_vec = []
    p_vec = []
    for a in BASE_ATOMS[n_type]:
        assert a in NORMALIZED_BASE[n_type]
        if a in p1:
            norm_vec.append(NORMALIZED_BASE[n_type][a])
            p_vec.append(p1[a])
    if len(p_vec)<3:
        return (None,None)
    sup = SVDSuperimposer()
    sup.set(np.array(norm_vec,'f'), np.array(p_vec,'f'))
    sup.run()
    (rot,tran) = sup.get_rotran()
    new_points = []
    for p in (p1,p2):
        atoms = list(p.keys())
        vec = []
        for a in atoms:
            vec.append(p[a])
        new_vec = np.dot(np.array(vec,'f'), rot)+tran
        new_points.append(dict(list(zip(atoms,new_vec))))
    return new_points

def _apply_rot_tran(points, rot, tran):
    all_atoms = list(points.keys())
    all_vec = np.array([points[x] for x in all_atoms])
    all_vec = np.dot(all_vec, rot)+tran
    return dict(list(zip(all_atoms,all_vec)))

def _superimpose_atoms(ref_points, points, atoms):
    if ref_points is None or points is None or atoms is None:
        return (None,None,None,None)
    ref_vec = []
    vec = []
    for a in atoms:
        if a in ref_points and a in points:
            ref_vec.append(ref_points[a])
            vec.append(points[a])
    if len(vec)<3:
        return (None,None,None,None)
    sup = SVDSuperimposer()
    sup.set(np.array(ref_vec,'f'), np.array(vec,'f'))
    sup.run()
    (rot,tran) = sup.get_rotran()
    rms = sup.get_rms()
    return (_apply_rot_tran(points,rot,tran), rot, tran, rms)

def fit_points(points, n_type):
    assert n_type[0] in NORMALIZED_BASE
    assert n_type[1] in NORMALIZED_BASE
    assert n_type[0] in BASE_ATOMS
    assert n_type[1] in BASE_ATOMS
    
    p1,p2 = points
    p1,rot,tran,_rms = _superimpose_atoms(NORMALIZED_BASE[n_type[0]], p1, BASE_ATOMS[n_type[0]])
    p2 = _apply_rot_tran(p2,rot,tran)
    
    res_p1 = dict([(k,np.array(v,'f')) for k,v in list(NORMALIZED_BASE[n_type[0]].items()) ])
    res_p2,_new_rot,_new_tran,_rms = _superimpose_atoms(p2, NORMALIZED_BASE[n_type[1]], BASE_ATOMS[n_type[1]])
    return (res_p1,res_p2)

def _rmsd_formula(points, ref_points, rot, tran, rmsd_atoms):
    (c1,c2) = points
    (p1,p2) = ref_points  
    (r1,r2) = rmsd_atoms
    vec1 = []
    vec2 = []
    for a in r1:
        if a not in p1 or a not in c1:
            return 1000.0
        vec1.append(p1[a])
        vec2.append(c1[a])
    for a in r2:
        if a not in p2 or a not in c2:
            return 1000.0
        vec1.append(p2[a])
        vec2.append(c2[a])
    vec1 = np.array(vec1,'f')
    vec2 = np.dot(np.array(vec2,'f'), rot)+tran
    diff = vec2-vec1
    return np.sqrt(sum(sum(diff*diff))/len(diff))

def rmsd_distance(points, ref_points, sup_atoms, rmsd_atoms=None, multiple_rmsd_variants=False):
    (c1,c2) = points
    (p1,p2) = ref_points  
    for atoms_list,c_res,p_res in (sup_atoms[0],c1,p1), (sup_atoms[1],c2,p2):
        for a in atoms_list:
            if a not in c_res or a not in p_res:
                return 1000.0
    ref_p = [p1[a] for a in sup_atoms[0]] + [p2[a] for a in sup_atoms[1]]
    cur_p = [c1[a] for a in sup_atoms[0]] + [c2[a] for a in sup_atoms[1]]
    sup = SVDSuperimposer()
    sup.set(np.array(ref_p,'f'), np.array(cur_p,'f'))
    sup.run()
    if rmsd_atoms is not None:
        (rot,tran) = sup.get_rotran()
        if multiple_rmsd_variants:
            return min([_rmsd_formula(points,ref_points,rot,tran,r) for r in rmsd_atoms])
        else:
            return _rmsd_formula(points,ref_points,rot,tran,rmsd_atoms)
    else:
        return sup.get_rms()

def bp_distance(points, ref_points, n_type, already_normalized=False, zfactor=0.25):
    return bp_distance3(points, ref_points, n_type, already_normalized, zfactor)

def bp_distance1(points, ref_points):
    return rmsd_distance(points, ref_points, (['C2','C4','C6'],[]), ([],['C2','C4','C6']))

def bp_distance2(points, ref_points, n_type, already_normalized=False, zfactor=0.25):
    rmsd_atoms = ["C1'",'C2','C4','C6']
    n_type0 = n_type[0].upper()
    if not already_normalized:
        _tmp1,n_points = normalize_points(points,n_type0)
        _tmp1,n_ref_points = normalize_points(ref_points,n_type0)
    else:
        n_points = points[1]
        n_ref_points = ref_points[1]
    if any([x not in n_points for x in rmsd_atoms]):
        return 999.0
    if any([x not in n_ref_points for x in rmsd_atoms]):
        return 999.0
    vec1 = np.array([n_points[x] for x in rmsd_atoms],'f')
    vec2 = np.array([n_ref_points[x] for x in rmsd_atoms],'f')
    vec1 = vec1 * (1.0,1.0,zfactor)
    vec2 = vec2 * (1.0,1.0,zfactor)
    diff = vec2-vec1
    return np.sqrt(sum(sum(diff*diff))/len(diff))

def bp_distance3(points, ref_points, n_type, already_normalized=False, zfactor=0.25):
    rmsd_atoms = ["C1'",'C2','C4','C6']
    n_type0 = n_type[0].upper()
    n_type1 = n_type[1].upper()
    if not already_normalized:
        _tmp1,n_points = normalize_points(points,n_type0)
        _tmp1,n_ref_points = normalize_points(ref_points,n_type0)
    else:
        n_points = points[1]
        n_ref_points = ref_points[1]
    if "C1'" not in n_points or "C1'" not in n_ref_points:
        rmsd_atoms = ['C2','C4','C6']
    if any([x not in n_points for x in rmsd_atoms]):
        return 999.0
    if any([x not in n_ref_points for x in rmsd_atoms]):
        return 999.0
    c1_vec = center_vector(n_points,n_type1)
    c2_vec = center_vector(n_ref_points,n_type1)
    if c1_vec is None or c2_vec is None:
        return 999.0

    vec1 = np.array([n_points[x] for x in rmsd_atoms],'f')
    vec2 = np.array([n_ref_points[x] for x in rmsd_atoms],'f')
    vec1 = (vec1-c1_vec) * (1.0,1.0,zfactor)
    vec2 = (vec2-c2_vec) * (1.0,1.0,zfactor)
    diff = vec2-vec1
    return 1.0*vector_length(c1_vec-c2_vec)+5.0*np.sqrt(sum(sum(diff*diff))/len(diff))


def bph_distance(points, ref_points, n_type):
    n_type0 = n_type[0].upper()
    _tmp1,n_points = normalize_points(points,n_type0)
    _tmp2,n_ref_points = normalize_points(ref_points,n_type0)
    if points is None or n_ref_points is None:
        return 999.0
    c = 0
    dist_sum = 0.0
    oxygens = list(set(PH_OXYGENS).intersection(set(n_points.keys())).intersection(set(n_ref_points)))
    if len(oxygens)<3 or 'P' not in n_points or 'P' not in n_ref_points:
        return 999.0
    res = 999.0
    p_diff = n_points['P']-n_ref_points['P']
    p_sum = sum(p_diff*p_diff)
    for p_oxygens in itertools.permutations(oxygens):
        curr_diff = np.array([n_points[o1]-n_ref_points[o2] for o1,o2 in zip(oxygens,p_oxygens)],'f')
        curr_sum = sum(sum(curr_diff*curr_diff))
        res = min(res, curr_sum+p_sum)
    return np.sqrt(res/(1+len(oxygens)))

def range_value_distance(range, value):
    x0 = range['avg_d']
    x = value
    diff = abs(x0-x)
    if x>=x0:
        delta = max(0.1, range['max_d']-x0)
    else:
        delta = max(0.1, x0-range['min_d'])
    return math.exp(-max(1,(diff*diff)/(delta*delta))+1)


def range_value_distance_old(ref_d, curr_d):
    curr_d = curr_distances[k1][k2]
    std_dev = max(0.1, ref_d['std_dev'])
    ref_d = ref_d['avg_d']
    diff = curr_d-ref_d
    return math.exp(-max(1,(diff*diff)/(std_dev*std_dev))+1)

def vector_length(v):
    if v is None:
        return None
    return np.sqrt(np.add.reduce(v*v))

# code snatched from Scientific.Geometry
def vector_angle(vec_a, vec_b):
    if vec_a is None or vec_b is None:
        return None
    cosa = np.add.reduce(vec_a*vec_b) / \
        np.sqrt(np.add.reduce(vec_a*vec_a) * \
        np.add.reduce(vec_b*vec_b))
    cosa = max(-1., min(1., cosa))
    return np.arccos(cosa) * ARCPI

def torsion_angle(a1,a2,a3,a4):
    # see http://www.iucr.org/__data/iucr/cif/software/pymmlib/pymmlib-1.0.0/mmLib/AtomMath.py
    # and https://github.com/pycogent/pycogent/blob/master/cogent/struct/dihedral.py
    a12 = np.array(a2,'f')-np.array(a1,'f')
    a23 = np.array(a3,'f')-np.array(a2,'f')
    a34 = np.array(a4,'f')-np.array(a3,'f')

    n12 = normal_vector(a12,a23)
    n34 = normal_vector(a23,a34)
    
    cross_n12_n34  = np.cross(n12, n34)
    direction      = np.dot(cross_n12_n34, a23)
    scalar_product = max(-1.0, min(1.0, np.dot(n12, n34)))
    
    angle = np.arccos(scalar_product) * ARCPI
    
    
    if direction<0.0:
        angle = -angle
        
    return angle

def center_vector(atoms,n_type):
    if atoms is None:
        return None
    normal_set = NORMAL_SUPPORT.get(n_type.upper())
    if normal_set is None:
        return None
    asum = np.array([0.0, 0.0, 0.0])
    for atomname in normal_set:
        if atomname not in atoms:
            return None
        asum += atoms[atomname]
    return asum / 6.0

def ribose_center_vector(atoms):
    if atoms is None:
        return None
    normal_set = ["C4'","C3'","C2'","C1'","O4'"]
    asum = np.array([0.0, 0.0, 0.0])
    for atomname in normal_set:
        if atomname not in atoms:
            return None
        asum += atoms[atomname]
    return asum / 5.0


def sugar_vector(atoms):
    # "O2'","C2'
    normal_set = ["O2'","O3'","O4'"]
    asum = np.array([0.0, 0.0, 0.0])
    for atomname in normal_set:
        if atomname not in atoms:
            return None
        asum += atoms[atomname]
    return asum / len(normal_set)

def phosphate_vector(atoms):
    if 'P' not in atoms:
        return None
    return atoms['P']

def gl_start_vector(atoms,n_type):
    if atoms is None:
        return None
    if n_type in ['C','U']:
        a1 = 'N1'
    elif n_type in ['A','G']:
        a1 = 'N9'
    else:
        return None
    if a1 not in atoms:
        return None
    return np.array(atoms[a1],'f')

def gl_vector(atoms,n_type):
    if n_type in ['C','U']:
        a1 = 'N1'
    elif n_type in ['A','G']:
        a1 = 'N9'
    else:
        return None
    a2 = "C1'"
    if a1 not in atoms or a2 not in atoms:
        return None
    return np.array(atoms[a2],'f')-np.array(atoms[a1],'f')

def normal_vector(v1,v2):
    normal = np.cross(v1, v2)
    normal = normal/vector_length(normal)
    return normal

# modified from Moderna code
def base_normal_vector(atoms, n_type):
    normal_set = NORMAL_SUPPORT.get(n_type.upper())
    if normal_set is None:
        return None
    for i in range(4):
        if normal_set[i] not in atoms:
            return None
    atoma = np.array(atoms[normal_set[1]],'f') - np.array(atoms[normal_set[0]],'f')
    atomb = np.array(atoms[normal_set[3]],'f') - np.array(atoms[normal_set[2]],'f')
    normal = np.cross(atoma, atomb)
    normal = normal/vector_length(normal)
    return normal

def residue_conformation(atoms, n_type, strict_definition=True):
    # see:
    # http://pl.scribd.com/doc/52213745/16/Conformations-About-the-Glycosidic-Bond
    # https://github.com/BGSU-RNA/FR3D/blob/master/FR3DSource/mSynList.m
    if n_type in ['C','U']:
        base_atoms = ['N1','C2']
    elif n_type in ['A','G']:
        base_atoms = ['N9','C4']
    else:
        return None
    r_atoms = ["C1'","O4'"]
    
    if any([a not in atoms for a in base_atoms+r_atoms]):
        return None

    chi = torsion_angle(atoms[r_atoms[1]], atoms[r_atoms[0]], atoms[base_atoms[0]], atoms[base_atoms[1]])
    if chi>=0.0 and chi<=90.0:
        return 1
    elif chi>=-90.0 and chi<0.0:
        if strict_definition:
            return -1
        else:
            return 1
    else: 
        return -1 # anti 

def _stacking_overlap(points, n_type, already_fitted=False):
    def _myDet(p, q, r):
        sum1 = q[0]*r[1] + p[0]*q[1] + r[0]*p[1]
        sum2 = q[0]*p[1] + r[0]*q[1] + p[0]*r[1]
    
        return sum1 - sum2
    
    def _isRightTurn(xxx_todo_changeme):
        (p, q, r) = xxx_todo_changeme
        assert p != q and q != r and p != r
        if _myDet(p, q, r) < 0:
            return 1
        else:
            return 0
    
    def _isPointInPolygon(r, P):
        # list of points is lister counter-closewise
        for i in range(len(P[:-1])):
            p, q = P[i], P[i+1]
            if _isRightTurn((p, q, r)):
                return False
        return True

    if not re.match('^[ACGU]{2}$',n_type):
        return 0

    p1,p2 = points
    if already_fitted:
        pp1,pp2 = p1,p2
    else:
        pp1,pp2 = fit_points((p1,p2), n_type)
        
    convex_points = {
        'A': list(zip(
                [ -2.100463,  0.000000,  4.271447,  1.920945,  0.230436, -2.100463],
                [  0.447145, -1.009320,  1.317924,  5.150733,  4.699718,  0.447145]
             )),
        'C': list(zip(
                [ -2.082733,  0.000000,  2.269450,  1.203833, -0.527970, -2.036772, -2.082733],
                [  0.123632, -1.010259, -0.120783,  4.411996,  4.602202,  2.647095,  0.123632]
             )),
        'G': list(zip(
                [ -2.101572,  0.000000,  4.872516,  5.295175,  3.613335,  1.396986, -0.751391, -2.101572],
                [  0.463584, -1.009529,  0.097781,  1.782283,  3.374547,  4.394213,  2.120132,  0.463584]
             )),
        'U': list(zip(
                [ -2.082780,  0.000000,  2.292490,  2.092152,  0.177156, -2.124577, -2.082780],
                [  0.111836, -1.008947, -0.048394,  2.445179,  4.020060,  2.616537,  0.111836]
             )),
    }
    res = []
    for atom_name,(x,y,z) in list(pp2.items()):
        if _isPointInPolygon((x,y), convex_points[n_type[0]]):
            res.append(atom_name)
    return len(res)

def stacking_overlap(points, n_type, already_fitted=False):
    (p1,p2) = points
    rev_n_type = n_type[::-1]
    return min(_stacking_overlap((p1,p2), n_type, already_fitted), _stacking_overlap((p2,p1), rev_n_type, already_fitted))

def base_min_dist(points, n_type, already_fitted=False):
    p1,p2 = points
    if already_fitted:
        pp1,pp2 = p1,p2
    else:
        pp1,pp2 = fit_points((p1,p2), n_type)
    return np.min( scipy.spatial.distance.cdist(np.array(list(pp1.values()),'f'), np.array(list(pp2.values()),'f')) )

def rotation_to_axis_angle(matrix):
    """Convert the rotation matrix into the axis-angle notation.
    source: http://svn.gna.org/svn/relax/tags/1.3.4/maths_fns/rotation_matrix.py

    Conversion equations
    ====================

    From Wikipedia (http://en.wikipedia.org/wiki/Rotation_matrix), the conversion is given by::

        x = Qzy-Qyz
        y = Qxz-Qzx
        z = Qyx-Qxy
        r = hypot(x,hypot(y,z))
        t = Qxx+Qyy+Qzz
        theta = atan2(r,t-1)

    @param matrix:  The 3x3 rotation matrix to update.
    @type matrix:   3x3 numpy array
    @return:    The 3D rotation axis and angle (in degrees).
    @rtype:     numpy 3D rank-1 array, float
    """

    # Axes.
    axis = np.zeros(3, np.float64)
    axis[0] = matrix[2,1] - matrix[1,2]
    axis[1] = matrix[0,2] - matrix[2,0]
    axis[2] = matrix[1,0] - matrix[0,1]

    # Angle.
    r = np.hypot(axis[0], np.hypot(axis[1], axis[2]))
    # inny wzor: r = math.sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2])
    t = matrix[0,0] + matrix[1,1] + matrix[2,2]
    theta = math.atan2(r, t-1)*ARCPI

    # Normalize the axis.
    axis = axis / r
        
    # Return the data.
    return axis, theta

def fr3d_axis_angle(rot):
    if rot[2][2]>0.0:
        axis, angle = rotation_to_axis_angle(rot)
    else:
        axis, angle = rotation_to_axis_angle(np.dot(rot, np.diag([-1,1,-1])))
    
    if axis[2] < 0: # based on FR3D rules
        axis = np.dot(axis,-1);
        angle = -angle

    if angle<-90.0:
        angle += 360
    return axis, angle

def doublet_params_bp_dict(points,n_type,type='bp'):
    """returns: 
      dist     -- distance between rings centers
      nn_ang   -- angle between normal vectors
      nn_ang_norm -- normalized angle between normal vectors = min(nn_ang,180-nn_ang)
      n1cc_ang -- angle between normal1 and centers vector
      n2cc_ang -- angle between normal2 and centers vector
      n12cc_ang - min(min(n1cc_ang,180.0-n1cc_ang), min(n2cc_ang,180.0-n2cc_ang))
      o_ang     - orientation of glycosidic bond
    """
    (p1,p2) = points
    c1_vec = center_vector(p1,n_type[0])
    n1_vec = base_normal_vector(p1,n_type[0])
    c2_vec = center_vector(p2,n_type[1])
    n2_vec = base_normal_vector(p2,n_type[1])
    if c1_vec is None or c2_vec is None:
        return None
    (fit1,fit2) = fit_points(points,n_type)
    if fit1 is None or fit2 is None:
        return None
    displ = gl_start_vector(fit2,n_type[1])-gl_start_vector(fit1,n_type[0])
    fit_n2_vec = base_normal_vector(fit2,n_type[1])
    min_dist = base_min_dist((fit1,fit2), n_type, already_fitted=True)
    
    # transform p1 to standard base
    _tmp, rot1, tran1,_rms = _superimpose_atoms(NORMALIZED_BASE.get(n_type[0]), p1, BASE_ATOMS.get(n_type[0]))
    # transform p2 to standard base
    _tmp, rot2, tran2,_rms = _superimpose_atoms(NORMALIZED_BASE.get(n_type[1]), p2, BASE_ATOMS.get(n_type[1]))
    # compute rotation matrix
    rot = np.dot(np.transpose(rot1), rot2)
    # compute axis and angle for rotation
    _axis, rot_ang = fr3d_axis_angle(rot)

    c_vec = c2_vec-c1_vec
    dist = vector_length(c_vec)
    nn_ang = vector_angle(n1_vec, n2_vec)
    nn_ang_norm = min(nn_ang,180.0-nn_ang)
    n1cc_ang = vector_angle(n1_vec, c_vec)
    n2cc_ang = vector_angle(n2_vec, c_vec)
    if n1cc_ang is not None and n2cc_ang is not None:
        n12cc_ang = min(min(n1cc_ang,180.0-n1cc_ang), min(n2cc_ang,180.0-n2cc_ang))
    else:
        n12cc_ang = None
    gl1_vec = gl_vector(p1,n_type[0])
    gl2_vec = gl_vector(p2,n_type[1])
    if gl1_vec is not None and gl2_vec is not None:
        o_ang = vector_angle(gl1_vec, gl2_vec)
        
        normal = np.cross(c_vec, n1_vec)
        v1 = np.dot(gl1_vec, normal)
        v2 = np.dot(gl2_vec, normal)
        if v1*v2 >= 0:
            orient = 1 # cis
        else:
            orient = -1 # tran
    else:
        o_ang = None
        orient = None
    strand_orient = None
    stack_orient = None
    dist_z = None
    if n1_vec is not None and n2_vec is not None:
        (pp1,pp2) = normalize_points(points,n_type[0])
        dist_z = abs((center_vector(pp2,n_type[1])-center_vector(pp1,n_type[0]))[2])
        nn1_vec = base_normal_vector(pp1,n_type[0])
        nn2_vec = base_normal_vector(pp2,n_type[1])
        v = nn1_vec[2]*nn2_vec[2]
        if v>=0:
            strand_orient = 1
        else:
            strand_orient = -1
    
    conf1 = residue_conformation(p1,n_type[0])
    conf2 = residue_conformation(p2,n_type[1])
    conf = None
    if conf1 is not None and conf2 is not None:
        if conf1==-1 and conf2==-1:
            conf = 0
        elif conf1==-1 and conf2==1:
            conf = 1
        elif conf1==1 and conf2==-1:
            conf = 2
        elif conf1==1 and conf2==1:
            conf = 3
    strand_orient_norm = strand_orient
    if conf is not None and strand_orient is not None:
        if conf in (1,2):
            strand_orient_norm = -strand_orient

    
    res = {'dist':dist,'dist_z':dist_z,
            'nn_ang':nn_ang,'nn_ang_norm':nn_ang_norm,
            'n1cc_ang':n1cc_ang,'n2cc_ang':n2cc_ang,'n12cc_ang':n12cc_ang,'o_ang': o_ang,
            'n2_z': fit_n2_vec[2],
            'rot_ang': rot_ang,
            'orient': orient,
            'strand_orient':strand_orient, 'strand_orient_norm': strand_orient_norm,
            'min_dist': min_dist,
            'conf':conf}
    if type in ['stacking','all','bp']: # TODO, rozdzielic type=bp od type=stacking (stacking_overlap) liczy sie dosyc dlugo!
        if nn_ang is not None and n1_vec is not None and n2_vec is not None:
            n1c2 = c_vec - n1_vec
            n1c2dist = vector_length(n1c2)
            is_up = n1c2dist < dist
    
            if nn_ang<90 and is_up==True:
                stack_orient = 0 # >>
            elif nn_ang<90 and is_up==False:
                stack_orient = 1 # <<
            elif nn_ang>=90 and is_up==True:
                stack_orient = 2 # ><
            elif nn_ang>=90 and is_up==False:
                stack_orient = 3 # <>

        s_overlap = None
        s_norm = None
        if n1_vec is not None and n2_vec is not None:
            s_overlap = stacking_overlap((fit1,fit2), n_type, already_fitted=True)
            s_norm = abs(fit_n2_vec[2])
        res['stack_orient']=stack_orient
        res['stack_overlap']=s_overlap
        res['stack_min_dist']=min_dist # TODO: usunac po 2013-02-18
        res['stack_norm']=s_norm
    return res

def doublet_params_bp_ph(points,n_type):
    (p1,p2) = points
    c1_vec = center_vector(p1,n_type[0])
    n1_vec = base_normal_vector(p1,n_type[0])
    ph2_vec = phosphate_vector(p2)
    if c1_vec is None or ph2_vec is None:
        return (None,None)
    ph_vec = ph2_vec-c1_vec
    dist_ph = vector_length(ph_vec)
    ph_ang = vector_angle(n1_vec, ph_vec)
    return (dist_ph, ph_ang)

def _compute_oxygen_interactions(p1,p2,ph_set,ox_set):
    res = {}
    i_oxygens = set()
    for m_atom,h_atom in ph_set:
        if m_atom in p1 and h_atom in p1:
            min_dist = None
            i_count = 0
            i_atom = None
            i_ang = None
            m_vec = np.array(p1[m_atom],'f')
            h_vec = np.array(p1[h_atom],'f')
            # print h_atom,h_vec
            for o_atom in ox_set:
                if o_atom not in p2:
                    continue
                o_vec = np.array(p2[o_atom],'f')
                cur_dist = vector_length(o_vec-m_vec)
                cur_ang = vector_angle(o_vec-h_vec, m_vec-h_vec)
                # print m_atom,h_atom,o_atom,cur_dist,cur_ang
                if min_dist is None or cur_dist<min_dist:
                    min_dist = cur_dist
                    i_atom = o_atom
                    i_ang = cur_ang
                if m_atom[0]=='C' and cur_dist<=4.0 and cur_ang>130:
                    i_count += 1
                    # print m_atom,h_atom,o_atom,cur_dist,cur_ang
                    i_oxygens.add(o_atom)
                    res['ii_%s_%s_%s'%(m_atom,h_atom,o_atom)]=1
                elif m_atom[0]=='N' and cur_dist<=3.5 and cur_ang>130:
                    i_count += 1
                    # print m_atom,h_atom,o_atom,cur_dist,cur_ang
                    i_oxygens.add(o_atom)
                    res['ii_%s_%s_%s'%(m_atom,h_atom,o_atom)]=1
            if min_dist is not None:
                res['dist_%s'%m_atom]=min_dist
                res['i_%s_%s'%(m_atom,h_atom)]=min(1,i_count)
                # res['i_atom_%s'%m_atom]=i_atom
    res['oxygens'] = list(i_oxygens)
    res['oxygens_count'] = len(i_oxygens)
    return res


def doublet_params_bp_br_dict(points,n_type):
    n_type0 = n_type[0].upper()
    br_set = BR_INTERACTIONS.get(n_type0)
    if br_set is None:
        return None
    
    p1 = NORMALIZED_BASE[n_type0]
    (_p1,p2) = normalize_points(points,n_type0)
    if p2 is None:
        return None
    
    c1_vec = center_vector(p1,n_type[0])
    n1_vec = base_normal_vector(p1,n_type[0])
    ph2_vec = phosphate_vector(p2)
    br_vec = ribose_center_vector(p2)

    if c1_vec is None or ph2_vec is None or br_vec is None:
        return None

    ph_vec = ph2_vec-c1_vec
    ph_ang = vector_angle(n1_vec, ph_vec)
    ph_dist = vector_length(ph_vec)
    ph_h = abs(ph_dist * np.cos(ph_ang/ARCPI))
    
    n1br_ang = vector_angle(n1_vec, br_vec-c1_vec)
    n1br_ang_norm = min(n1br_ang, 180-n1br_ang)
    
    res = {'ph_h': ph_h,'ph_dist':ph_dist, 'br_dist': vector_length(br_vec-c1_vec), 'n1br_ang': n1br_ang, 'n1br_ang_norm': n1br_ang_norm}
    res.update(_compute_oxygen_interactions(p1,p2,br_set,BR_OXYGENS))
    return res

def doublet_params_bp_ph_dict(points,n_type):
    n_type0 = n_type[0].upper()
    ph_set = PH_INTERACTIONS.get(n_type0)
    if ph_set is None:
        return None
    p1 = NORMALIZED_BASE[n_type0]
    (_p1,p2) = normalize_points(points,n_type0)
    if p2 is None:
        return None

    c1_vec = center_vector(p1,n_type[0])
    n1_vec = base_normal_vector(p1,n_type[0])
    ph2_vec = phosphate_vector(p2)

    if c1_vec is None or ph2_vec is None:
        return None

    ph_vec = ph2_vec-c1_vec
    ph_ang = vector_angle(n1_vec, ph_vec)
    ph_dist = vector_length(ph_vec)
    ph_h = abs(ph_dist * np.cos(ph_ang/ARCPI))
    
    res = {'ph_h': ph_h,'ph_dist': ph_dist, 'n1ph_ang':ph_ang}
    res.update(_compute_oxygen_interactions(p1,p2,ph_set,PH_OXYGENS))
    return res

def doublet_params_dict(points,n_type,type="bp"):
    if type=='bp' or type=='stacking':
        return doublet_params_bp_dict(points, n_type,type=type)
    elif type=='base-phosphate':
        return doublet_params_bp_ph_dict(points,n_type)
    elif type=='base-ribose':
        return doublet_params_bp_br_dict(points,n_type)
    else:
        raise Exception("unknown param type: %s" % type)
        
def expected_strand_orient(desc):
    if desc in ['WW_cis','HH_cis','SS_cis','WS_cis','SW_cis']:
        return -1 # anti-parallel
    elif desc in ['WW_tran','HH_tran','SS_tran','WS_tran','SW_tran']:
        return 1 # parallel
    elif desc in ['WH_cis','HW_cis','HS_cis','SH_cis']:
        return 1 # parallel
    elif desc in ['WH_tran','HW_tran','HS_tran','SH_tran']:
        return -1 # anti-parallel
    else:
        raise Exception("Unknown desc: %s" % desc)
        
def expected_stack_orient(desc):
    if desc=='>>':
        return 0
    elif desc=='<<':
        return 1
    elif desc=='><':
        return 2
    elif desc=='<>':
        return 3
    else:
        raise Exception("Unknown desc: %s" % desc)


######################################################################################
######################################################################################
######################################################################################


def method_cache(f):
    from functools import wraps
    # print("cacher called")
    cache_attr = "_"+f.__name__
    @wraps(f)
    def wrapped(self,*args, **kwds):
        # print "wrapped called %s, args=%s"%(f.func_name,args)
        if not hasattr(self,cache_attr):
            # print("calculating and caching result")
            setattr(self,cache_attr,f(self,*args, **kwds))
        else:
            # print "using cache"
            pass
        return getattr(self,cache_attr)
    return wrapped

class RotatedResidue:
    def __init__(self,r,rot,tran):
        self.r = r
        self.rot = rot
        self.tran = tran
        
    @property
    def n_type(self):
        return self.r.n_type

    @property
    @method_cache
    def points(self):
        return _apply_rot_tran(self.r.points_dict, self.rot, self.tran)
        
    @property
    def center(self):
        v = self.r.center
        if v is None:
            return None
        return np.dot(v,self.rot)+self.tran

    @property
    def normal_vector(self):
        return np.dot(self.r.normal_vector,self.rot) # uwaga! tutaj bez translacji

    @property
    def gl_vector(self):
        return np.dot(self.r.gl_vector,self.rot)+self.tran

    @property
    def ph_vec(self):
        v = self.r.ph_vec
        if v is None:
            return None
        return np.dot(v,self.rot)+self.tran

    @property
    def br_vec(self):
        v = self.r.br_vec
        if v is None:
            return None
        return np.dot(v,self.rot)+self.tran

    @property
    def conf(self):
        return self.r.conf

class Residue:

    def __init__(self,r_id,n_type,points,repair_missing=True):
        self.r_id = r_id
        tmp = r_id.split(":")
        if len(tmp)==1:
            self.pdb_id,self.chain,self.res_num = (None,tmp[0][0],tmp[0][1:])
        else:
            self.pdb_id,self.chain,self.res_num = (tmp[0],tmp[1][0],tmp[1][1:])
        self.n_type = n_type
        self.points_dict = dict([(k,np.array(p,'f')) for k,p in list(points.items())])
        
        norm2points, self.fit_rot2, self.fit_tran2, rms = _superimpose_atoms(self.points_dict, NORMALIZED_BASE.get(n_type), BASE_ATOMS.get(n_type))

        common_base_atoms = set(BASE_ATOMS.get(self.n_type,[])).intersection(list(self.points_dict.keys()))
        if repair_missing and norm2points is not None:
            if rms>0.3:
                for k,v in list(norm2points.items()):
                    self.points_dict[k] = v
            elif norm2points is not None and len(common_base_atoms)!=len(BASE_ATOMS.get(self.n_type,[])):
                for k,v in list(norm2points.items()):
                    if k not in self.points_dict:
                        self.points_dict[k] = v
        self.normalized_points, self.fit_rot, self.fit_tran, _rms = _superimpose_atoms(NORMALIZED_BASE.get(n_type), self.points_dict, BASE_ATOMS.get(n_type))
        
        
        #self.fit_points2 = NORMALIZED_BASE[n_type]
    
    @property
    def points(self):
        return self.points_dict
    
    @property
    def fit_points(self):
        return NORMALIZED_RESIDUES[self.n_type].points
    
    @property
    def fit(self):
        return NORMALIZED_RESIDUES[self.n_type]
    
    @property
    @method_cache
    def center(self):
        return center_vector(self.points, self.n_type)

    @property
    @method_cache
    def normal_vector(self):
        return base_normal_vector(self.points, self.n_type)
        
    @property
    @method_cache
    def gl_vector(self):
        return gl_vector(self.points, self.n_type)

    @property
    @method_cache
    def conf(self):
        return residue_conformation(self.points, self.n_type)

    @property
    @method_cache
    def ph_vec(self):
        return phosphate_vector(self.points)
    
    @property
    @method_cache
    def br_vec(self):
        return ribose_center_vector(self.points)

NORMALIZED_RESIDUES = dict([(n_type,Residue("A1",n_type,p)) for n_type,p in list(NORMALIZED_BASE.items())])
    
class Doublet:

    def __init__(self,d_id,res1,res2):
        self.d_id = d_id
        tmp = d_id.split(":")
        if len(tmp)==2:
            self.pdb_id,self.res1_id,self.res2_id = (None,tmp[0],tmp[1])
        else:
            self.pdb_id,self.res1_id,self.res2_id = (tmp[0],tmp[1],tmp[2])
        self.res1 = res1
        self.res2 = res2
        self.n_type = res1.n_type+res2.n_type
        
        if self.res2.fit_rot2 is not None and self.res1.fit_rot is not None:
            rot = np.dot(self.res2.fit_rot2, self.res1.fit_rot)
            tran = np.dot(self.res2.fit_tran2, self.res1.fit_rot)+self.res1.fit_tran

            self.fit1 = self.res1.fit
            self.fit2 = RotatedResidue(self.res2.fit, rot, tran) # self.res1.fit_rot2, self.res1.fit_tran2)

            self.normalized1 = RotatedResidue(self.res1, self.res1.fit_rot, self.res1.fit_tran)
            self.normalized2 = RotatedResidue(self.res2, self.res1.fit_rot, self.res1.fit_tran)
        else:
            self.fit1 = None
            self.fit2 = None
            self.normalized1 = None
            self.normalized2 = None
    
    def get_all_params(self):
        fields = (
                    "nn_ang","nn_ang_norm","n1cc_ang","n2cc_ang","n12cc_ang","n1ph_ang",
                    "o_ang","orient","stack_orient","dist","min_dist","dist_z","strand_orient",
                    "conf","strand_orient_norm","stack_norm",
                    "n2_z","rot_ang","stack_overlap",
                    "ph_dist","ph_ang","ph_h","ph_info",
                    "br_dist","n1br_ang_norm","br_info",
                 )
        return dict([(k,getattr(self,k)) for k in fields])
    
    @property
    def center_vector(self):
        """vector between centers of r1 & r2"""
        if not hasattr(self,'_center_vector'):
            if self.res1.center is None or self.res2.center is None:
                self._center_vector = None
            else:
                self._center_vector = self.res2.center - self.res1.center
        return self._center_vector

    @property
    def nn_ang(self):
        if not hasattr(self,'_nn_ang'):
            self._nn_ang = vector_angle(self.res1.normal_vector, self.res2.normal_vector)
        return self._nn_ang
    
    @property
    def nn_ang_norm(self):
        if not hasattr(self,'_nn_ang_norm'):
            if self.nn_ang is not None:
                self._nn_ang_norm = min(self.nn_ang, 180.0-self.nn_ang)
            else:
                self._nn_ang_norm = None
        return self._nn_ang_norm
    
    @property
    @method_cache
    def n1cc_ang(self):
        return vector_angle(self.res1.normal_vector, self.center_vector)

    @property
    @method_cache
    def n1ph_ang(self):
        return vector_angle(self.fit1.normal_vector, self.ph_vec)

    @property
    @method_cache
    def n2cc_ang(self):
        return vector_angle(self.res2.normal_vector, self.center_vector)
    
    @property
    @method_cache
    def n12cc_ang(self):
        if self.n1cc_ang is None or self.n2cc_ang is None:
            return None
        n1cc = self.n1cc_ang
        n2cc = self.n2cc_ang
        return min(min(n1cc,180.0-n1cc), min(n2cc,180.0-n2cc))

    @property
    @method_cache
    def o_ang(self):
        """angle between glycosidic vectors"""
        return vector_angle(self.res1.gl_vector, self.res2.gl_vector)

    @property
    def orient(self):
        """kind of orientation, TODO: do we still use this?"""
        if not hasattr(self,'_orient'):
            if self.center_vector is not None and self.res1.normal_vector is not None \
               and self.res1.gl_vector is not None and self.res2.gl_vector is not None:
                normal = np.cross(self.center_vector, self.res1.normal_vector)
                v1 = np.dot(self.res1.gl_vector, normal)
                v2 = np.dot(self.res2.gl_vector, normal)
                if v1*v2 >= 0:
                    self._orient = 1 # cis
                else:
                    self._orient = -1 # tran
            else:
                self._orient = None
        return self._orient

    @property
    @method_cache
    def stack_orient(self):
        if self.center_vector is None or self.res1.normal_vector is None \
           or self.dist is None or self.nn_ang is None:
           return None
        n1c2 = self.center_vector - self.res1.normal_vector
        n1c2dist = vector_length(n1c2)
        is_up = n1c2dist < self.dist
        if self.nn_ang<90 and is_up==True:
            return 0 # >>
        elif self.nn_ang<90 and is_up==False:
            return 1 # <<
        elif self.nn_ang>=90 and is_up==True:
            return 2 # ><
        elif self.nn_ang>=90 and is_up==False:
            return 3 # <>
    
    @property
    def dist(self):
        if not hasattr(self,'_dist'):
            self._dist = vector_length(self.center_vector)
        return self._dist
    
    @property
    @method_cache
    def min_dist(self):
        """minimal distance between atoms in r1 & r2"""
        p1 = np.array(list(self.fit1.points.values()),'f')
        p2 = np.array(list(self.fit2.points.values()),'f')
        return np.min( scipy.spatial.distance.cdist(p1, p2) )
        
    @property
    @method_cache
    def dist_z(self):
        c1 = self.normalized1.center
        c2 = self.normalized2.center
        if c1 is None or c2 is None:
            return None
        return abs((c1-c2)[2])
    
    @property
    @method_cache
    def strand_orient(self):
        nn1_vec = self.normalized1.normal_vector
        nn2_vec = self.normalized2.normal_vector
        if nn1_vec is None or nn2_vec is None:
            return None
        v = nn1_vec[2]*nn2_vec[2]
        if v>=0:
            return 1
        else:
            return -1
    
    @property
    def conf(self):
        conf1 = self.res1.conf
        conf2 = self.res2.conf
        if conf1 is None or conf2 is None:
            return None
        
        if conf1==-1 and conf2==-1:
            return 0
        elif conf1==-1 and conf2==1:
            return 1
        elif conf1==1 and conf2==-1:
            return 2
        elif conf1==1 and conf2==1:
            return 3
        else:
            return None
            
    @property
    @method_cache
    def strand_orient_norm(self):
        res = self.strand_orient
        if res is not None and self.conf in (1,2):
            res = -res
        return res
        
    @property
    @method_cache
    def stack_norm(self):
        if self.fit_n2_vec is not None:
            return abs(self.fit_n2_vec[2])
        else:
            return None
    
    @property
    @method_cache
    def fit_n2_vec(self):
        return base_normal_vector(self.fit2.points,self.res2.n_type)

    @property
    @method_cache
    def n2_z(self):
        if self.fit_n2_vec is not None:
            return self.fit_n2_vec[2]
        else:
            return None

    @property
    @method_cache
    def rot_ang(self):
        rot1 = self.res1.fit_rot
        rot2 = self.res2.fit_rot
        if rot1 is None or rot2 is None:
            return None
        
        rot = np.dot(np.transpose(rot1), rot2)
        # compute axis and angle for rotation
        _axis, rot_ang = fr3d_axis_angle(rot)
        return rot_ang

    @property
    @method_cache
    def stack_overlap(self):
        res = stacking_overlap((self.fit1.points,self.fit2.points), self.n_type, already_fitted=True)
        return res

    # base-phosphate

    @property
    @method_cache
    def ph_vec(self):
        if self.fit1 is None or self.normalized2 is None:
            return None
        center1 = self.fit1.center
        ph_vec2 = self.normalized2.ph_vec
        if ph_vec2 is None or center1 is None:
            return None
        res = ph_vec2-center1
        return res

    @property
    @method_cache
    def ph_dist(self):
        return vector_length(self.ph_vec)

    @property
    @method_cache
    def ph_ang(self):
        if self.fit1 is None:
            return None
        normal1 = self.fit1.normal_vector
        ph_vec = self.ph_vec
        if ph_vec is None or normal1 is None:
            return None
        return vector_angle(normal1, ph_vec)
    
    @property
    @method_cache
    def ph_h(self):
        ph_ang = self.ph_ang
        ph_dist = self.ph_dist
        if ph_ang is None or ph_dist is None:
            return None
        ph_h = abs(ph_dist * np.cos(ph_ang/ARCPI))
        return ph_h

    @property
    @method_cache
    def ph_info(self):
        if self.fit1 is None or self.normalized2 is None:
            return None
        ph_set = PH_INTERACTIONS.get(self.res1.n_type)
        return _compute_oxygen_interactions(self.fit1.points,self.normalized2.points,ph_set,PH_OXYGENS)
        
    # base-ribose

    @property
    @method_cache
    def br_vec(self):
        if self.fit1 is None or self.normalized2 is None:
            return None
        center1 = self.fit1.center
        br_vec2 = self.normalized2.br_vec
        if br_vec2 is None or center1 is None:
            return None
        return br_vec2-center1

    @property
    @method_cache
    def br_dist(self):
        br_vec = self.br_vec
        if br_vec is None:
            return None
        return vector_length(br_vec)
    
    @property
    @method_cache
    def n1br_ang(self):
        normal1 = self.fit1.normal_vector
        br_vec = self.br_vec
        return vector_angle(normal1, br_vec)

    @property
    @method_cache
    def n1br_ang_norm(self):
        n1br_ang = self.n1br_ang
        if n1br_ang is None:
            return None
        return min(n1br_ang, 180-n1br_ang)
        
    @property
    @method_cache
    def br_info(self):
        if self.fit1 is None or self.normalized2 is None:
            return None
        br_set = BR_INTERACTIONS.get(self.res1.n_type)
        return _compute_oxygen_interactions(self.fit1.points,self.normalized2.points,br_set,BR_OXYGENS)
    
    
    # OUTDATED
    @property
    def stack_min_dist(self):
        return self.min_dist

