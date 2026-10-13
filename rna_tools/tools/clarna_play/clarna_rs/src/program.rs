//! ClaRNA classifier programs (ClaRNAlib/clarna.py: ContactProgram, MyExpr*), loaded from classifier.*.json.gz.
//!
//! Options as used by rna_clarna_run.py: disable_align_to and disable_compare_distances are on,
//! so ALIGN_TO* only sets rmsd = 1000 and COMPARE_DISTANCES only sets dist_score = 1000.

use crate::geom::{dist, f32r, V3};
use crate::residue::{center_vector, Doublet, Points};
use serde_json::Value as J;

#[derive(Clone, Debug, PartialEq)]
pub enum Val {
    None,
    Num(f64),
    Str(String),
    Bool(bool),
    V(V3),
}

impl Val {
    fn truthy(&self) -> bool {
        match self {
            Val::None => false,
            Val::Num(x) => *x != 0.0,
            Val::Str(s) => !s.is_empty(),
            Val::Bool(b) => *b,
            Val::V(_) => true,
        }
    }
    fn num(&self) -> Option<f64> {
        match self {
            Val::Num(x) => Some(*x),
            Val::Bool(b) => Some(if *b { 1.0 } else { 0.0 }),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Attr {
    NType,
    SimpleNType,
    Rmsd,
    BpDistance,
    DistanceToPoints,
    DistScore,
    Param(&'static str),
    Norm2Center,
}

#[derive(Debug)]
pub enum Arg {
    Attr(Attr),
    Lit(Val),
    Expr(Box<Expr>),
}

/// reference points of COMPUTE_BP_DISTANCE, prepared once at load time
#[derive(Debug)]
pub struct BpRef {
    center: Option<V3>,
    with_c1: bool,
    atoms: Vec<Option<V3>>, // C1', C2, C4, C6 (float32, as np.array(..., 'f'))
}

#[derive(Debug)]
pub enum Expr {
    If(Vec<Expr>, Vec<Expr>),
    Align(Option<String>),
    ComputeBpDistance(BpRef, Option<String>),
    CompareDistances(Option<String>),
    DistanceToPoints(Arg, Vec<V3>),
    ComputeDistanceToPoints(Arg, Vec<V3>),
    Le(Arg, Arg),
    Lt(Arg, Arg),
    Eq(Arg, Arg),
    Or(Vec<Expr>),
    InRange(Arg, Arg, Arg),
    Return(String, Option<Arg>),
    ReturnIf(Arg, f64, f64, String),
}

#[derive(Clone, Debug)]
pub struct MyResult {
    pub desc: String,
    pub score: f64,
    pub real_score: f64,
    pub reference_id: Option<String>,
}

pub struct Program {
    pub commands: Vec<Expr>,
}

fn attr_of(name: &str) -> Attr {
    match name {
        "n_type" => Attr::NType,
        "simple_n_type" => Attr::SimpleNType,
        "rmsd" => Attr::Rmsd,
        "bp_distance" => Attr::BpDistance,
        "distance_to_points" => Attr::DistanceToPoints,
        "dist_score" => Attr::DistScore,
        "norm2_center" => Attr::Norm2Center,
        _ => {
            if name.starts_with("d_param_br_") || name.starts_with("d_param_ph_") || name.starts_with("norm1_") || name.starts_with("norm2_") {
                panic!("attribute not implemented: {}", name);
            }
            if let Some(k) = name.strip_prefix("d_param_") {
                const KNOWN: [&str; 13] = [
                    "dist", "nn_ang", "nn_ang_norm", "n1cc_ang", "n2cc_ang", "n12cc_ang", "min_dist", "strand_orient",
                    "stack_norm", "rot_ang", "stack_overlap", "stack_orient", "stack_min_dist",
                ];
                let k = if k == "stack_min_dist" { "min_dist" } else { k };
                match KNOWN.iter().find(|x| **x == k) {
                    Some(x) => Attr::Param(x),
                    None => panic!("doublet parameter not implemented: {}", k),
                }
            } else {
                panic!("attribute not implemented: {}", name)
            }
        }
    }
}

fn lit(j: &J) -> Val {
    match j {
        J::Null => Val::None,
        J::Bool(b) => Val::Bool(*b),
        J::Number(n) => Val::Num(n.as_f64().unwrap()),
        J::String(s) => Val::Str(s.clone()),
        _ => panic!("unsupported literal {}", j),
    }
}

/// MyExpression.arg_value: "$attr" -> attribute, list -> expression (only LE loads these), else literal
fn arg(j: &J, allow_expr: bool) -> Arg {
    match j {
        J::String(s) if s.len() > 1 && s.starts_with('$') => Arg::Attr(attr_of(&s[1..])),
        J::Array(_) if allow_expr => Arg::Expr(Box::new(load_expr(j))),
        _ => Arg::Lit(lit(j)),
    }
}

fn points(j: &J) -> Vec<V3> {
    j.as_array()
        .unwrap()
        .iter()
        .map(|p| {
            let p = p.as_array().unwrap();
            [p[0].as_f64().unwrap(), p[1].as_f64().unwrap(), p[2].as_f64().unwrap()]
        })
        .collect()
}

fn opt_str(j: Option<&J>) -> Option<String> {
    j.and_then(|x| x.as_str()).map(|s| s.to_string())
}

fn bp_ref(j: &J) -> BpRef {
    // ref_points[1]: atom -> [x, y, z]
    let d = j.as_array().unwrap()[1].as_object().unwrap();
    let p = Points(
        d.iter()
            .map(|(k, v)| {
                let v = v.as_array().unwrap();
                (k.clone(), [v[0].as_f64().unwrap(), v[1].as_f64().unwrap(), v[2].as_f64().unwrap()])
            })
            .collect(),
    );
    BpRef {
        // the ring atoms are the same set for all bases, so the center does not depend on n_type
        center: center_vector(&p, 'A'),
        with_c1: p.get("C1'").is_some(),
        atoms: ["C1'", "C2", "C4", "C6"].iter().map(|a| p.get(a).map(f32r)).collect(),
    }
}

pub fn load_expr(j: &J) -> Expr {
    let p = j.as_array().unwrap();
    let op = p[0].as_str().unwrap();
    match op {
        "IF" => Expr::If(load_program(&p[1]), load_program(&p[2])),
        "ALIGN_TO" => Expr::Align(opt_str(p.get(3))),
        "ALIGN_TO2" => Expr::Align(opt_str(p.get(4))),
        "ALIGN_TO_WITH_RMSD_ATOMS" | "ALIGN_TO_WITH_MULTIPLE_RMSD_ATOMS" => Expr::Align(if p.len() == 6 { opt_str(p.get(4)) } else { None }),
        "COMPUTE_BP_DISTANCE" => Expr::ComputeBpDistance(bp_ref(&p[1]), opt_str(p.get(2))),
        "COMPARE_DISTANCES" => Expr::CompareDistances(opt_str(p.get(2))),
        "DISTANCE_TO_POINTS" => Expr::DistanceToPoints(arg(&p[1], false), points(&p[2])),
        "COMPUTE_DISTANCE_TO_POINTS" => Expr::ComputeDistanceToPoints(arg(&p[1], false), points(&p[2])),
        "LE" => Expr::Le(arg(&p[1], true), arg(&p[2], true)),
        "LT" => Expr::Lt(arg(&p[1], false), arg(&p[2], false)),
        "GE" => Expr::Le(arg(&p[2], false), arg(&p[1], false)),
        "GT" => Expr::Lt(arg(&p[2], false), arg(&p[1], false)),
        "EQ" => Expr::Eq(arg(&p[1], false), arg(&p[2], false)),
        "OR" => Expr::Or(p[1..].iter().map(load_expr).collect()),
        "IN_RANGE" => Expr::InRange(arg(&p[1], false), arg(&p[2], false), arg(&p[3], false)),
        "RETURN" => Expr::Return(p[1].as_str().unwrap().to_string(), p.get(2).map(|x| arg(x, false))),
        "RETURN_IF" => Expr::ReturnIf(
            arg(&p[1], false),
            p[2].as_f64().unwrap(),
            p[3].as_f64().unwrap(),
            p[4].as_str().unwrap().to_string(),
        ),
        _ => panic!("unknown expression {}", op),
    }
}

fn load_program(j: &J) -> Vec<Expr> {
    j.as_array().unwrap().iter().map(load_expr).collect()
}

impl Program {
    pub fn load(path: &str) -> Program {
        let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{}: {}", path, e));
        let mut s = String::new();
        use std::io::Read;
        flate2::read::GzDecoder::new(f).read_to_string(&mut s).unwrap();
        let j: J = serde_json::from_str(&s).unwrap();
        Program { commands: load_program(&j) }
    }

    /// ContactProgram.value: all results of the program for a doublet (None if there are none)
    pub fn value(&self, d: &Doublet) -> Option<Vec<MyResult>> {
        let mut st = State { d, rmsd: None, bp_distance: None, distance_to_points: None, dist_score: None, reference: None };
        let mut res: Option<Vec<MyResult>> = None;
        for c in &self.commands {
            collect(&mut res, st.eval(c));
        }
        res
    }
}

enum Out {
    Val(Val),
    Res(Option<Vec<MyResult>>),
    One(MyResult),
}

fn collect(res: &mut Option<Vec<MyResult>>, o: Out) {
    match o {
        Out::Res(Some(v)) => res.get_or_insert_with(Vec::new).extend(v),
        Out::One(r) => res.get_or_insert_with(Vec::new).push(r),
        _ => {}
    }
}

fn truthy(o: &Out) -> bool {
    match o {
        Out::Val(v) => v.truthy(),
        Out::Res(r) => r.is_some(),
        Out::One(_) => true,
    }
}

/// ContactProgramState (with a Doublet)
struct State<'a, 'b> {
    d: &'b Doublet<'a>,
    rmsd: Option<f64>,
    bp_distance: Option<f64>,
    distance_to_points: Option<f64>,
    dist_score: Option<f64>,
    reference: Option<String>,
}

fn opt(x: Option<f64>) -> Val {
    x.map(Val::Num).unwrap_or(Val::None)
}

impl<'a, 'b> State<'a, 'b> {
    fn attr(&self, a: Attr) -> Val {
        match a {
            Attr::NType => Val::Str(self.d.n_type.clone()),
            Attr::SimpleNType => {
                let s: Vec<&str> = self
                    .d
                    .n_type
                    .chars()
                    .filter_map(|c| match c {
                        'A' | 'G' => Some("Pu"),
                        'C' | 'U' => Some("Py"),
                        _ => None,
                    })
                    .collect();
                Val::Str(s.join("-"))
            }
            Attr::Rmsd => opt(self.rmsd),
            Attr::BpDistance => opt(self.bp_distance),
            Attr::DistanceToPoints => opt(self.distance_to_points),
            Attr::DistScore => opt(self.dist_score),
            Attr::Param(k) => opt(self.d.param(k)),
            Attr::Norm2Center => self.d.norm2_center.map(Val::V).unwrap_or(Val::None),
        }
    }

    fn arg(&mut self, a: &Arg) -> Val {
        match a {
            Arg::Attr(x) => self.attr(*x),
            Arg::Lit(v) => v.clone(),
            Arg::Expr(e) => match self.eval(e) {
                Out::Val(v) => v,
                Out::Res(None) => Val::None,
                _ => Val::Bool(true),
            },
        }
    }

    fn nearest(&mut self, a: &Arg, pts: &[V3]) -> f64 {
        match self.arg(a) {
            Val::V(p) => pts.iter().map(|q| dist(p, *q)).fold(f64::INFINITY, f64::min),
            _ => 999.0,
        }
    }

    /// distances.bp_distance3 with already_normalized=True
    fn bp_distance(&self, r: &BpRef) -> f64 {
        let n_points = match &self.d.norm2_points {
            Some(p) => p,
            None => return 999.0, // python would fail here (normalized2 is None)
        };
        let nt1 = self.d.n_type.chars().nth(1).unwrap_or('?').to_ascii_uppercase();
        let idx: &[usize] = if n_points.get("C1'").is_some() && r.with_c1 { &[0, 1, 2, 3] } else { &[1, 2, 3] };
        const NAMES: [&str; 4] = ["C1'", "C2", "C4", "C6"];
        let mut v1 = vec![];
        let mut v2 = vec![];
        for &i in idx {
            match (n_points.get(NAMES[i]), r.atoms[i]) {
                (Some(a), Some(b)) => {
                    v1.push(f32r(a));
                    v2.push(b);
                }
                _ => return 999.0,
            }
        }
        let c1 = match center_vector(n_points, nt1) {
            Some(c) => c,
            None => return 999.0,
        };
        let c2 = match (crate::residue::normal_support(nt1), r.center) {
            (Some(_), Some(c)) => c,
            _ => return 999.0,
        };
        let mut s = 0.0;
        for (a, b) in v1.iter().zip(&v2) {
            let mut d = [0.0; 3];
            for k in 0..3 {
                let z = if k == 2 { 0.25 } else { 1.0 };
                d[k] = (b[k] - c2[k]) * z - (a[k] - c1[k]) * z;
            }
            s += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        }
        dist(c1, c2) + 5.0 * (s / v1.len() as f64).sqrt()
    }

    fn eval(&mut self, e: &Expr) -> Out {
        match e {
            Expr::If(conds, then) => {
                for c in conds {
                    let o = self.eval(c);
                    if !truthy(&o) {
                        return Out::Res(None);
                    }
                }
                let mut res = None;
                for c in then {
                    let o = self.eval(c);
                    collect(&mut res, o);
                }
                Out::Res(res)
            }
            Expr::Or(conds) => {
                for c in conds {
                    let o = self.eval(c);
                    if truthy(&o) {
                        return Out::Val(Val::Bool(true));
                    }
                }
                Out::Val(Val::Bool(false))
            }
            Expr::Align(refid) => {
                self.rmsd = Some(1000.0);
                self.reference = refid.clone();
                Out::Res(None)
            }
            Expr::ComputeBpDistance(r, refid) => {
                self.reference = refid.clone();
                self.bp_distance = Some(self.bp_distance(r));
                Out::Res(None)
            }
            Expr::CompareDistances(refid) => {
                self.dist_score = Some(1000.0);
                self.reference = refid.clone();
                Out::Res(None)
            }
            Expr::DistanceToPoints(a, pts) => Out::Val(Val::Num(self.nearest(a, pts))),
            Expr::ComputeDistanceToPoints(a, pts) => {
                self.distance_to_points = Some(self.nearest(a, pts));
                Out::Res(None)
            }
            Expr::Le(a, b) => {
                let (x, y) = (self.arg(a), self.arg(b));
                Out::Val(Val::Bool(match (x.num(), y.num()) {
                    (Some(x), Some(y)) => x <= y,
                    _ => false,
                }))
            }
            Expr::Lt(a, b) => {
                let (x, y) = (self.arg(a), self.arg(b));
                Out::Val(Val::Bool(match (x.num(), y.num()) {
                    (Some(x), Some(y)) => x < y,
                    _ => false,
                }))
            }
            Expr::Eq(a, b) => {
                let (x, y) = (self.arg(a), self.arg(b));
                let eq = match (x.num(), y.num()) {
                    (Some(x), Some(y)) => x == y,
                    _ => x == y,
                };
                Out::Val(Val::Bool(eq))
            }
            Expr::InRange(a, lo, hi) => {
                let (x, l, h) = (self.arg(a), self.arg(lo), self.arg(hi));
                Out::Val(Val::Bool(match (x.num(), l.num(), h.num()) {
                    (Some(x), Some(l), Some(h)) => x >= l && x <= h,
                    _ => false,
                }))
            }
            Expr::Return(desc, score) => {
                let s = score.as_ref().and_then(|a| self.arg(a).num());
                Out::One(MyResult {
                    desc: desc.clone(),
                    score: s.unwrap_or(1.0),
                    real_score: s.unwrap_or(1000.0),
                    reference_id: self.reference.clone(),
                })
            }
            Expr::ReturnIf(test, ok, fuzzy, desc) => {
                let s = match self.arg(test).num() {
                    Some(s) => s,
                    None => return Out::Res(None),
                };
                if s <= *ok {
                    let r = (1.0 - 0.5 * s / ok.max(0.00000001)).clamp(0.0, 1.0);
                    return Out::One(MyResult { desc: desc.clone(), score: r, real_score: s, reference_id: self.reference.clone() });
                }
                if s <= *fuzzy {
                    let bd = (fuzzy - ok).max(0.1);
                    return Out::One(MyResult {
                        desc: format!("?{}", desc),
                        score: 0.5 - 0.5 * (s - ok) / bd,
                        real_score: s,
                        reference_id: self.reference.clone(),
                    });
                }
                Out::Res(None)
            }
        }
    }
}
