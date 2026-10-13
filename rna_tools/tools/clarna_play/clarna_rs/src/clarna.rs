//! The ClaRNA pipeline (rna_clarna_run.py: extract_close_doublets_from_fobject, start_clarna,
//! clarna._find_contacts_using_library, SeeClarna.make_structure_dict / display_bbstack).

use crate::geom::{dist, V3};
use crate::pdb;
use crate::program::{MyResult, Program};
use crate::residue::{Doublet, NormalizedResidues, Points, Residue};

pub const MAX_RES_DIST: f64 = 7.0;

/// an entry of `s`: [chain.id, resname, resnum, Residue, id]
pub struct SEntry {
    pub resname: String,
    pub id: String,
    pub res: Residue,
}

pub struct Structure {
    pub s: Vec<SEntry>,
    pub close_doublets: Vec<(usize, usize)>,
}

/// rna_clarna_run.extract_close_doublets_from_fobject (+ utils.compute_close_doublets, neighbours sorted)
pub fn extract(models: &[pdb::Model]) -> Structure {
    let mut s: Vec<SEntry> = vec![];
    let mut close = vec![];
    for model in models {
        // residues of this model (None: not A/C/G/U), indices are model-local (as in ClaRNA)
        let mut residues: Vec<Option<Vec<(String, V3)>>> = vec![];
        for chain in &model.chains {
            for r in chain.residues() {
                let resname = r.resname.to_uppercase();
                // simplify_residue: duplicate atom ids -> None (residue skipped)
                let mut seen = std::collections::HashSet::new();
                if !r.atoms.iter().all(|a| seen.insert(a.id.clone())) {
                    continue;
                }
                let atoms: Vec<(String, V3)> = r.atoms.iter().map(|a| (a.id.clone(), a.coord())).collect();
                let id = format!("{}{}", chain.id, r.resseq);
                let rn = resname.trim().to_string();
                residues.push(if ["A", "C", "G", "U"].contains(&rn.as_str()) { Some(atoms.clone()) } else { None });
                s.push(SEntry { resname: rn.clone(), id, res: Residue::new(&rn, Points(atoms)) });
            }
        }
        for (i, j) in compute_close_doublets(&residues, MAX_RES_DIST) {
            close.push((i, j));
            close.push((j, i));
        }
    }
    Structure { s, close_doublets: close }
}

fn compute_close_doublets(residues: &[Option<Vec<(String, V3)>>], max_distance: f64) -> Vec<(usize, usize)> {
    let mut num = vec![];
    let mut atoms: Vec<&Vec<(String, V3)>> = vec![];
    let mut refs: Vec<V3> = vec![];
    for (i, r) in residues.iter().enumerate() {
        let r = match r {
            Some(r) => r,
            None => continue,
        };
        // a_dict: last value wins for duplicate names (not possible here)
        let get = |k: &str| r.iter().rev().find(|(n, _)| n == k).map(|(_, v)| *v);
        if let Some(p) = ["C1'", "C2", "C4", "C6", "P"].iter().find_map(|a| get(a)) {
            num.push(i);
            atoms.push(r);
            refs.push(p);
        }
    }
    let mut result = vec![];
    for ii in 0..num.len() {
        let i = num[ii];
        for jj in 0..num.len() {
            if dist(refs[ii], refs[jj]) > max_distance + 15.0 {
                continue;
            }
            let j = num[jj];
            // `j >= len(residues_atoms_np)` compares a residue index with a count: ClaRNA quirk, kept
            if j <= i || j >= atoms.len() || residues[j].is_none() {
                continue;
            }
            let mut d = f64::INFINITY;
            for (_, a) in atoms[ii] {
                for (_, b) in atoms[jj] {
                    d = d.min(dist(*a, *b));
                }
            }
            if d >= 0.0 && d <= max_distance {
                result.push((i, j));
            }
        }
    }
    result
}

#[derive(Clone, Debug)]
pub struct Contact {
    pub n_type: String,
    pub desc: String,
    pub full_desc: String,
    pub weight: f64,
    pub reference: Option<String>,
    pub is_alt: bool,
}

#[derive(Clone, Debug)]
pub enum Edge {
    Dist { n_type: String },
    Contact(Contact),
}

/// networkx MultiDiGraph with insertion-ordered nodes, neighbours and edge keys
#[derive(Default)]
pub struct Graph {
    nodes: Vec<String>,
    adj: Vec<Vec<(usize, Vec<(String, Edge)>)>>,
}

impl Graph {
    fn node(&mut self, n: &str) -> usize {
        match self.nodes.iter().position(|x| x == n) {
            Some(i) => i,
            None => {
                self.nodes.push(n.to_string());
                self.adj.push(vec![]);
                self.nodes.len() - 1
            }
        }
    }
    fn add_edge(&mut self, u: &str, v: &str, key: Option<&str>, e: Edge) {
        let (u, v) = (self.node(u), self.node(v));
        let nb = &mut self.adj[u];
        let k = match nb.iter().position(|(x, _)| *x == v) {
            Some(k) => k,
            None => {
                nb.push((v, vec![]));
                nb.len() - 1
            }
        };
        let keys = &mut nb[k].1;
        let key = match key {
            Some(k) => k.to_string(),
            None => {
                // networkx: the lowest unused integer key
                let mut n = keys.len();
                while keys.iter().any(|(x, _)| *x == n.to_string()) {
                    n += 1;
                }
                n.to_string()
            }
        };
        match keys.iter_mut().find(|(x, _)| *x == key) {
            Some(x) => x.1 = e,
            None => keys.push((key, e)),
        }
    }
    pub fn edges(&self) -> impl Iterator<Item = (&str, &str, &Edge)> {
        self.adj.iter().enumerate().flat_map(move |(u, nb)| {
            nb.iter().flat_map(move |(v, keys)| keys.iter().map(move |(_, e)| (self.nodes[u].as_str(), self.nodes[*v].as_str(), e)))
        })
    }
}

/// DoubletDescTool.get_desc_category (non-detailed)
fn category(desc: &str) -> &'static str {
    if desc.is_empty() {
        return "unclassified";
    }
    let d = if desc.len() > 1 && (desc.starts_with('?') || desc.starts_with('!')) { &desc[1..] } else { desc };
    let b = d.as_bytes();
    let n = b.iter().take_while(|c| b"sSWH".contains(c)).take(2).count();
    if n >= 1 && (d[n..].starts_with("_cis") || d[n..].starts_with("_tran")) {
        return "bp";
    }
    if b.len() >= 2 && b"<>".contains(&b[0]) && b"<>".contains(&b[1]) {
        return "stacking";
    }
    let n = b.iter().take_while(|c| b"SWH".contains(c)).take(2).count();
    if n >= 1 && d[n..].starts_with('_') {
        let rest = &d[n + 1..];
        let nd = rest.bytes().take_while(|c| c.is_ascii_digit()).count();
        if nd >= 1 && rest[nd..].starts_with("BP") {
            return "base-phosphate";
        }
        if nd >= 1 && rest[nd..].starts_with("BR") {
            return "base-ribose";
        }
    }
    "other"
}

/// clarna.aggregate_results
fn aggregate(res: &[MyResult]) -> Vec<MyResult> {
    let mut a: Vec<(String, MyResult)> = vec![];
    for r in res {
        if r.desc.starts_with('!') {
            continue;
        }
        let d = r.desc.replace('?', "");
        match a.iter_mut().find(|(k, _)| *k == d) {
            Some(e) => {
                if r.score > e.1.score {
                    e.1 = r.clone();
                }
            }
            None => a.push((d, r.clone())),
        }
    }
    let mut v: Vec<MyResult> = a.into_iter().map(|(_, r)| r).collect();
    v.sort_by(|x, y| (-x.score).partial_cmp(&-y.score).unwrap().then(x.real_score.partial_cmp(&y.real_score).unwrap()));
    v
}

fn simple_desc(x: &str) -> String {
    x.replace('?', "").replace('!', "").split(' ').next().unwrap().to_string()
}

pub struct ClaRNA {
    libs: Vec<Program>,
    nres: NormalizedResidues,
}

impl ClaRNA {
    pub fn new(lib_paths: &[String]) -> ClaRNA {
        ClaRNA { libs: lib_paths.iter().map(|p| Program::load(p)).collect(), nres: NormalizedResidues::new() }
    }

    /// clarna._find_contacts_using_library
    pub fn classify(&self, st: &Structure) -> Graph {
        let mut g = Graph::default();
        for e in &st.s {
            g.node(&e.id);
        }
        for &(i, j) in &st.close_doublets {
            let (si, sj) = (&st.s[i], &st.s[j]);
            let n_type = format!("{}{}", si.resname, sj.resname);
            g.add_edge(&si.id, &sj.id, None, Edge::Dist { n_type: n_type.clone() });
            for lib in &self.libs {
                let d = Doublet::new(&si.res, &sj.res, &self.nres);
                let res = lib.value(&d);
                // divide_results_by_cat
                let mut cats: Vec<(&str, Vec<MyResult>)> =
                    vec![("bp", vec![]), ("stacking", vec![]), ("base-phosphate", vec![]), ("base-ribose", vec![])];
                for r in res.unwrap_or_default() {
                    let c = category(&r.desc);
                    match cats.iter_mut().find(|(k, _)| *k == c) {
                        Some(x) => x.1.push(r),
                        None => match cats.iter_mut().find(|(k, _)| *k == "unk") {
                            Some(x) => x.1.push(r),
                            None => cats.push(("unk", vec![r])),
                        },
                    }
                }
                let mut cur_dict: Vec<MyResult> = vec![];
                let mut alt_dict: Vec<MyResult> = vec![];
                for (_, results) in &cats {
                    if results.is_empty() {
                        continue;
                    }
                    let agg = aggregate(results);
                    let cur = match agg.first() {
                        Some(c) => c.clone(),
                        None => continue,
                    };
                    let d1 = simple_desc(&cur.desc);
                    // sorted(aggregate_results(results), key=-score): stable, agg is already in that order
                    let alt = agg.iter().find(|r| simple_desc(&r.desc) != d1);
                    if let Some(alt) = alt {
                        if cur.desc.contains('?') && alt.score > 0.7 * cur.score {
                            alt_dict.push(alt.clone());
                        }
                    }
                    cur_dict.push(cur);
                }
                for (is_alt, dict) in [(false, &cur_dict), (true, &alt_dict)] {
                    for r in dict.iter() {
                        let d = r.desc.split(' ').next().unwrap().to_string();
                        g.add_edge(
                            &si.id,
                            &sj.id,
                            Some(&d),
                            Edge::Contact(Contact {
                                n_type: n_type.clone(),
                                desc: d.clone(),
                                full_desc: r.desc.clone(),
                                weight: r.score,
                                reference: r.reference_id.clone(),
                                is_alt,
                            }),
                        );
                    }
                }
            }
        }
        g
    }
}

/// SeeClarna.make_structure_dict + display_bbstack, for edges with weight > min_score
pub fn display(g: &Graph, min_score: f64) -> String {
    // (uc, ux, vc, vx) -> ("bp": (un, vn), [(t, w)])
    let mut dict: Vec<((char, i64, char, i64), (char, char), Vec<(String, f64)>)> = vec![];
    for (u, v, e) in g.edges() {
        let c = match e {
            Edge::Contact(c) if c.weight > min_score => c,
            _ => continue,
        };
        let (xc, yc) = (u.chars().next().unwrap(), v.chars().next().unwrap());
        let (xr, yr): (i64, i64) = (u[xc.len_utf8()..].trim().parse().unwrap(), v[yc.len_utf8()..].trim().parse().unwrap());
        let nt: Vec<char> = c.n_type.chars().collect();
        let (xn, yn) = (nt[0], nt[1]);
        let mut t = c.desc.clone();
        let key;
        let bp;
        if yr < xr {
            key = (yc, yr, xc, xr);
            bp = (yn, xn);
            // reverse the orientation of the interaction
            if t.starts_with('<') || t.starts_with('>') {
                if t == "<<" {
                    t = ">>".to_string();
                } else if t == ">>" {
                    t = "<<".to_string();
                }
            } else {
                let ch: Vec<char> = t.chars().collect();
                t = format!("{}{}{}", ch[1], ch[0], ch[2..].iter().collect::<String>());
            }
        } else {
            key = (xc, xr, yc, yr);
            bp = (xn, yn);
        }
        let w = c.weight;
        match dict.iter_mut().find(|(k, _, _)| *k == key) {
            Some((_, _, types)) => match types.iter_mut().find(|(x, _)| *x == t) {
                Some(x) => x.1 = (x.1 + w) * 0.5,
                None => types.push((t, w)),
            },
            None => dict.push((key, bp, vec![(t, w)])),
        }
    }
    // bubbleSort on ux: a stable sort
    dict.sort_by_key(|(k, _, _)| k.1);
    let mut out = String::new();
    for ((uc, ux, vc, vx), (un, vn), types) in &dict {
        out.push_str(&format!("{} {:4}   {} {:4}", uc, ux, vc, vx));
        out.push_str(&format!("          bp {} {}", un, vn));
        for (t, w) in types {
            out.push_str(&format!("    {:>20} {:8.4}", t, w));
        }
        out.push('\n');
    }
    out
}
