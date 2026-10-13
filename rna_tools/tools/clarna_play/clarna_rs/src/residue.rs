//! ClaRNA's Residue and Doublet (ClaRNAlib/distances.py), only what the bp and stacking classifiers use.

use crate::geom::*;

/// atom name -> coordinates, in insertion order (like a Python dict)
#[derive(Clone, Debug, Default)]
pub struct Points(pub Vec<(String, V3)>);

impl Points {
    pub fn get(&self, k: &str) -> Option<V3> {
        self.0.iter().find(|(n, _)| n == k).map(|(_, v)| *v)
    }
    pub fn set(&mut self, k: &str, v: V3) {
        match self.0.iter_mut().find(|(n, _)| n == k) {
            Some(e) => e.1 = v,
            None => self.0.push((k.to_string(), v)),
        }
    }
    /// x·rot + tran; `round`: the Python code works on a float32 array here
    pub fn apply(&self, rot: &M3, tran: V3, round: bool) -> Points {
        Points(
            self.0
                .iter()
                .map(|(n, v)| {
                    let x = add(vm(*v, rot), tran);
                    (n.clone(), if round { f32r(x) } else { x })
                })
                .collect(),
        )
    }
}

pub fn normalized_base(n: char) -> Option<&'static [(&'static str, V3)]> {
    match n {
        'A' => Some(&NB_A),
        'C' => Some(&NB_C),
        'U' => Some(&NB_U),
        'G' => Some(&NB_G),
        _ => None,
    }
}

pub fn base_atoms(n: char) -> &'static [&'static str] {
    match n {
        'A' => &["N9", "C4", "N3", "N1", "C6", "N6", "C8", "C5", "C2", "N7"],
        'C' => &["N1", "C2", "O2", "N3", "C4", "N4", "C6", "C5"],
        'U' => &["N1", "C2", "O2", "N3", "C4", "O4", "C6", "C5"],
        'G' => &["N9", "C4", "N3", "N1", "C6", "O6", "C8", "C5", "C2", "N7", "N2"],
        _ => &[],
    }
}

pub fn normal_support(n: char) -> Option<[&'static str; 6]> {
    match n {
        'C' | 'U' => Some(["N1", "C2", "N3", "C4", "C5", "C6"]),
        'G' | 'A' => Some(["N1", "C2", "C4", "N3", "C5", "C6"]),
        _ => None,
    }
}

/// center of the base ring (distances.center_vector)
pub fn center_vector(p: &Points, n: char) -> Option<V3> {
    let set = normal_support(n)?;
    let mut s = [0.0; 3];
    for a in set {
        s = add(s, p.get(a)?);
    }
    Some(scale(s, 1.0 / 6.0))
}

/// distances.base_normal_vector
pub fn base_normal_vector(p: &Points, n: char) -> Option<V3> {
    let set = normal_support(n)?;
    let a0 = p.get(set[0])?;
    let a1 = p.get(set[1])?;
    let a2 = p.get(set[2])?;
    let a3 = p.get(set[3])?;
    let nv = cross(sub(a1, a0), sub(a3, a2));
    Some(scale(nv, 1.0 / length(nv)))
}

/// distances._superimpose_atoms(ref_points, points, atoms): moves `points` onto `ref_points`
fn superimpose(reference: &Points, points: &Points, atoms: &[&str]) -> Option<(M3, V3, f64)> {
    let mut rv = vec![];
    let mut pv = vec![];
    for a in atoms {
        if let (Some(r), Some(p)) = (reference.get(a), points.get(a)) {
            rv.push(r);
            pv.push(p);
        }
    }
    if pv.len() < 3 {
        return None;
    }
    Some(kabsch(&rv, &pv))
}

fn nb_points(n: char) -> Option<Points> {
    normalized_base(n).map(|v| Points(v.iter().map(|(k, p)| (k.to_string(), *p)).collect()))
}

pub struct Residue {
    pub n_type: String,
    pub nt: char, // n_type as a single letter ('?' if not A/C/G/U)
    pub points: Points,
    pub fit_rot: Option<(M3, V3)>,  // residue -> normalized frame
    pub fit_rot2: Option<(M3, V3)>, // normalized base -> residue frame
    pub center: Option<V3>,
    pub normal: Option<V3>,
}

impl Residue {
    pub fn new(n_type: &str, points: Points) -> Residue {
        let nt = if n_type.len() == 1 { n_type.chars().next().unwrap() } else { '?' };
        let mut points = Points(points.0.into_iter().map(|(k, v)| (k, f32r(v))).collect());
        let nb = nb_points(nt);
        let atoms = base_atoms(nt);
        let mut fit_rot2 = None;
        if let Some(nb) = &nb {
            if let Some((rot2, tran2, rms)) = superimpose(&points, nb, atoms) {
                fit_rot2 = Some((rot2, tran2));
                // the idealized base (python floats -> float64 array) placed onto this residue
                let norm2points = nb.apply(&rot2, tran2, false);
                let common = atoms.iter().filter(|a| points.get(a).is_some()).count();
                if rms > 0.3 {
                    for (k, v) in norm2points.0 {
                        points.set(&k, v);
                    }
                } else if common != atoms.len() {
                    for (k, v) in norm2points.0 {
                        if points.get(&k).is_none() {
                            points.set(&k, v);
                        }
                    }
                }
            }
        }
        let fit_rot = nb.as_ref().and_then(|nb| superimpose(nb, &points, atoms)).map(|(r, t, _)| (r, t));
        let center = center_vector(&points, nt);
        let normal = base_normal_vector(&points, nt);
        Residue { n_type: n_type.to_string(), nt, points, fit_rot, fit_rot2, center, normal }
    }
}

/// NORMALIZED_RESIDUES[n_type]: a Residue built from the idealized base
pub struct NormalizedResidues {
    res: Vec<(char, Residue)>,
}

impl NormalizedResidues {
    pub fn new() -> NormalizedResidues {
        let res = ['A', 'C', 'G', 'U']
            .iter()
            .map(|&n| {
                let p = nb_points(n).unwrap();
                (n, Residue::new(&n.to_string(), p))
            })
            .collect();
        NormalizedResidues { res }
    }
    pub fn get(&self, n: char) -> Option<&Residue> {
        self.res.iter().find(|(k, _)| *k == n).map(|(_, r)| r)
    }
}

/// a pair of residues with lazily computed parameters (distances.Doublet)
pub struct Doublet<'a> {
    pub res1: &'a Residue,
    pub res2: &'a Residue,
    pub n_type: String,
    fit1: Option<&'a Residue>,
    fit2_points: Option<Points>,
    norm_rot: Option<(M3, V3)>, // res1.fit_rot: maps both residues into res1's normalized frame
    pub norm2_points: Option<Points>,
    pub norm2_center: Option<V3>,
    cache: std::cell::RefCell<Vec<(&'static str, Option<f64>)>>,
}

impl<'a> Doublet<'a> {
    pub fn new(res1: &'a Residue, res2: &'a Residue, nres: &'a NormalizedResidues) -> Doublet<'a> {
        let n_type = format!("{}{}", res1.n_type, res2.n_type);
        let (mut fit1, mut fit2_points, mut norm_rot, mut norm2_points, mut norm2_center) = (None, None, None, None, None);
        if let (Some((r2, t2)), Some((r1, t1))) = (res2.fit_rot2, res1.fit_rot) {
            let rot = mm(&r2, &r1);
            let tran = add(vm(t2, &r1), t1);
            fit1 = nres.get(res1.nt);
            fit2_points = nres.get(res2.nt).map(|f| f.points.apply(&rot, tran, true));
            norm_rot = Some((r1, t1));
            norm2_points = Some(res2.points.apply(&r1, t1, true));
            norm2_center = res2.center.map(|c| add(vm(c, &r1), t1));
        }
        Doublet { res1, res2, n_type, fit1, fit2_points, norm_rot, norm2_points, norm2_center, cache: Default::default() }
    }

    fn center_vector(&self) -> Option<V3> {
        Some(sub(self.res2.center?, self.res1.center?))
    }

    fn nn_ang(&self) -> Option<f64> {
        Some(angle(self.res1.normal?, self.res2.normal?))
    }

    fn dist(&self) -> Option<f64> {
        self.center_vector().map(length)
    }

    fn n1cc_ang(&self) -> Option<f64> {
        Some(angle(self.res1.normal?, self.center_vector()?))
    }

    fn n2cc_ang(&self) -> Option<f64> {
        Some(angle(self.res2.normal?, self.center_vector()?))
    }

    /// a doublet parameter, as getattr(doublet, key)
    pub fn param(&self, key: &'static str) -> Option<f64> {
        if let Some((_, v)) = self.cache.borrow().iter().find(|(k, _)| *k == key) {
            return *v;
        }
        let v = self.compute(key);
        self.cache.borrow_mut().push((key, v));
        v
    }

    fn compute(&self, key: &str) -> Option<f64> {
        match key {
            "dist" => self.dist(),
            "nn_ang" => self.nn_ang(),
            "nn_ang_norm" => self.nn_ang().map(|a| a.min(180.0 - a)),
            "n1cc_ang" => self.n1cc_ang(),
            "n2cc_ang" => self.n2cc_ang(),
            "n12cc_ang" => {
                let (n1, n2) = (self.n1cc_ang()?, self.n2cc_ang()?);
                Some(n1.min(180.0 - n1).min(n2.min(180.0 - n2)))
            }
            "min_dist" => {
                let p1 = &self.fit1?.points;
                let p2 = self.fit2_points.as_ref()?;
                let mut m = f64::INFINITY;
                for (_, a) in &p1.0 {
                    for (_, b) in &p2.0 {
                        m = m.min(dist(*a, *b));
                    }
                }
                Some(m)
            }
            "strand_orient" => {
                let (r1, _) = self.norm_rot?;
                let n1 = vm(self.res1.normal?, &r1);
                let n2 = vm(self.res2.normal?, &r1);
                Some(if n1[2] * n2[2] >= 0.0 { 1.0 } else { -1.0 })
            }
            "stack_norm" => {
                let v = base_normal_vector(self.fit2_points.as_ref()?, self.res2.nt)?;
                Some(v[2].abs())
            }
            "rot_ang" => {
                let (rot1, _) = self.res1.fit_rot?;
                let (rot2, _) = self.res2.fit_rot?;
                Some(fr3d_axis_angle(&mm(&transpose(&rot1), &rot2)))
            }
            "stack_overlap" => {
                let p1 = &self.fit1?.points;
                let p2 = self.fit2_points.as_ref()?;
                let n: Vec<char> = self.n_type.chars().collect();
                if n.len() != 2 || !n.iter().all(|c| "ACGU".contains(*c)) {
                    return Some(0.0);
                }
                let a = stacking_overlap(p2, n[0]);
                let b = stacking_overlap(p1, n[1]);
                Some(a.min(b) as f64)
            }
            "stack_orient" => {
                let cv = self.center_vector()?;
                let n1 = self.res1.normal?;
                let d = self.dist()?;
                let nn = self.nn_ang()?;
                let is_up = length(sub(cv, n1)) < d;
                Some(match (nn < 90.0, is_up) {
                    (true, true) => 0.0,
                    (true, false) => 1.0,
                    (false, true) => 2.0,
                    (false, false) => 3.0,
                })
            }
            _ => panic!("doublet parameter not implemented: {}", key),
        }
    }
}

/// distances.fr3d_axis_angle -> angle
fn fr3d_axis_angle(rot: &M3) -> f64 {
    let m = if rot[2][2] > 0.0 {
        *rot
    } else {
        mm(rot, &[[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
    };
    let axis = [m[2][1] - m[1][2], m[0][2] - m[2][0], m[1][0] - m[0][1]];
    let r = axis[0].hypot(axis[1].hypot(axis[2]));
    let t = m[0][0] + m[1][1] + m[2][2];
    let mut angle = r.atan2(t - 1.0) * 180.0 / std::f64::consts::PI;
    // axis/r: for r == 0 numpy gives nan and the sign test is false
    if r != 0.0 && axis[2] / r < 0.0 {
        angle = -angle;
    }
    if angle < -90.0 {
        angle += 360.0;
    }
    angle
}

/// number of atoms of `points` whose (x, y) lies inside the convex hull of base `n` (distances._stacking_overlap)
fn stacking_overlap(points: &Points, n: char) -> usize {
    let poly: &[(f64, f64)] = match n {
        'A' => &[(-2.100463, 0.447145), (0.000000, -1.009320), (4.271447, 1.317924), (1.920945, 5.150733), (0.230436, 4.699718), (-2.100463, 0.447145)],
        'C' => &[(-2.082733, 0.123632), (0.000000, -1.010259), (2.269450, -0.120783), (1.203833, 4.411996), (-0.527970, 4.602202), (-2.036772, 2.647095), (-2.082733, 0.123632)],
        'G' => &[(-2.101572, 0.463584), (0.000000, -1.009529), (4.872516, 0.097781), (5.295175, 1.782283), (3.613335, 3.374547), (1.396986, 4.394213), (-0.751391, 2.120132), (-2.101572, 0.463584)],
        'U' => &[(-2.082780, 0.111836), (0.000000, -1.008947), (2.292490, -0.048394), (2.092152, 2.445179), (0.177156, 4.020060), (-2.124577, 2.616537), (-2.082780, 0.111836)],
        _ => return 0,
    };
    let det = |p: (f64, f64), q: (f64, f64), r: (f64, f64)| {
        let s1 = q.0 * r.1 + p.0 * q.1 + r.0 * p.1;
        let s2 = q.0 * p.1 + r.0 * q.1 + p.0 * r.1;
        s1 - s2
    };
    points
        .0
        .iter()
        .filter(|(_, v)| {
            let r = (v[0], v[1]);
            (0..poly.len() - 1).all(|i| det(poly[i], poly[i + 1], r) >= 0.0)
        })
        .count()
}

macro_rules! nb {
    ($($n:expr => [$x:expr, $y:expr, $z:expr]),* $(,)?) => { [$(($n, [$x, $y, $z])),*] };
}

// NORMALIZED_BASE (distances.py), in the same key order
static NB_A: [(&str, V3); 15] = nb![
    "N7" => [-0.7729614406105143, 2.1201797334590102, 0.0], "H9" => [0.0, -1.0093196646068083, 0.0],
    "1H6" => [1.9209447414524516, 5.150733145313508, 0.0], "H2" => [4.2714470975421825, 1.3179237109436448, 0.0],
    "H8" => [-2.1004626204274777, 0.44714489628299675, 0.0], "2H6" => [0.23043597510468228, 4.699718288132125, 0.0],
    "N9" => [0.0, 0.0, 0.0], "C8" => [-1.0899103005137063, 0.8335872143873861, 0.0],
    "C5" => [0.6063291460346252, 2.1057981896577656, 0.0], "C2" => [3.202584142514681, 1.51596183391717, 0.0],
    "N1" => [2.8648446154295764, 2.824390566297474, 0.0], "N3" => [2.409029202573013, 0.4396032568253272, 0.0],
    "C6" => [1.563644845822768, 3.1435063393974403, 0.0], "N6" => [1.2031042415497286, 4.447033338301733, 0.0],
    "C4" => [1.1184904351494842, 0.8031271160408657, 0.0],
];
static NB_C: [(&str, V3); 13] = nb![
    "1H4" => [-0.5279700289179345, 4.60220182612121, 0.0], "H1" => [0.0, -1.010258706209454, 0.0],
    "H6" => [-2.082732930782349, 0.12363197089251765, 0.0], "2H4" => [1.2038334533450203, 4.411995903837219, 0.0],
    "C6" => [-1.162053283211799, 0.6983551068687612, 0.0], "H5" => [-2.0367724443934563, 2.647095419696937, 0.0],
    "C2" => [1.28933438636058, 0.5992529922671946, 0.0], "N1" => [0.0, 0.0, 0.0],
    "N3" => [1.312136920098135, 1.977774853061991, 0.0], "N4" => [0.28135537444115655, 4.009202043572646, 0.0],
    "O2" => [2.2694499775195744, -0.12078289494958472, 0.0], "C5" => [-1.1289628161952552, 2.058236947549645, 0.0],
    "C4" => [0.18336964933177424, 2.654626399935204, 0.0],
];
static NB_U: [(&str, V3); 12] = nb![
    "N3" => [1.1953166017937347, 1.9738920601305123, 0.0], "H3" => [2.0921523976816854, 2.445178553234257, 0.0],
    "H1" => [0.0, -1.008947070985391, 0.0], "H6" => [-2.082780056506479, 0.11183564649907562, 0.0],
    "C6" => [-1.1777603199476516, 0.7096949325723031, 0.0], "H5" => [-2.1245767823028237, 2.616536750689694, 0.0],
    "C2" => [1.2602866366915804, 0.58925012217277, 0.0], "N1" => [0.0, 0.0, 0.0],
    "O4" => [0.17715643875791365, 4.020059648793736, 0.0], "O2" => [2.292490396543795, -0.04839432941443833, 0.0],
    "C5" => [-1.1961299637178628, 2.0625141391615487, 0.0], "C4" => [0.060003385436139034, 2.8083560986778724, 0.0],
];
static NB_G: [(&str, V3); 16] = nb![
    "1H2" => [4.872516232460849, 0.09778091003631939, 0.0], "N7" => [-0.7513911514240377, 2.1201318166099155, 0.0],
    "C2" => [3.2594124387538095, 1.333440824418896, 0.0], "H9" => [0.0, -1.0095288732671288, 0.0],
    "H1" => [3.613335050886427, 3.3745474050771, 0.0], "2H2" => [5.295175245866893, 1.7822826455998741, 0.0],
    "H8" => [-2.101572214696437, 0.4635842689564793, 0.0], "N9" => [0.0, 0.0, 0.0],
    "C8" => [-1.0886393160814454, 0.841883274438163, 0.0], "N1" => [2.8974071829603005, 2.6580453952999106, 0.0],
    "O6" => [1.3969862391333745, 4.394213054302783, 0.0], "N2" => [4.595966940484209, 1.0633818952456446, 0.0],
    "N3" => [2.4038183708450203, 0.34099632542447367, 0.0], "C6" => [1.5671222764357484, 3.191068234453134, 0.0],
    "C5" => [0.6252146444859401, 2.096901945838508, 0.0], "C4" => [1.117587071308519, 0.7907046767118864, 0.0],
];
