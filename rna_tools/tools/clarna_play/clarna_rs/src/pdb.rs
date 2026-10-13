//! PDB reader that mirrors Bio.PDB.PDBParser (permissive mode) closely enough for ClaRNA:
//! models, chains in order of first appearance, residue ids (hetero flag, resseq, icode),
//! altlocs (highest occupancy wins), duplicate atoms/residues, and coordinates stored as float32.

use crate::geom::{f32r, V3};

#[derive(Clone, Debug)]
pub struct Atom {
    pub id: String,
    // disordered atom: (altloc, coord, occupancy); `selected` is the visible one
    children: Vec<(char, V3, f64, String)>,
    selected: usize,
    last_occ: f64,
    disordered: bool,
}

impl Atom {
    pub fn coord(&self) -> V3 {
        self.children[self.selected].1
    }
    fn altloc(&self) -> char {
        self.children[self.selected].0
    }
}

#[derive(Clone, Debug)]
pub struct Residue {
    pub resname: String,
    pub resseq: i64,
    pub atoms: Vec<Atom>,
}

impl Residue {
    fn find(&self, name: &str) -> Option<usize> {
        self.atoms.iter().position(|a| a.id == name)
    }
}

/// a chain entry: a residue, or a disordered residue (several residues with different names)
#[derive(Clone, Debug)]
struct ResEntry {
    id: (String, i64, char),
    variants: Vec<Residue>,
    selected: usize,
}

#[derive(Clone, Debug)]
pub struct Chain {
    pub id: char,
    entries: Vec<ResEntry>,
}

impl Chain {
    /// residues as Bio iterates them (selected variant of disordered residues)
    pub fn residues(&self) -> impl Iterator<Item = &Residue> {
        self.entries.iter().map(|e| &e.variants[e.selected])
    }
}

#[derive(Clone, Debug, Default)]
pub struct Model {
    pub chains: Vec<Chain>,
}

/// where atoms of the current residue go
enum Cur {
    None,
    Entry(usize, usize), // chain index, entry index
    Orphan,              // a residue Bio created but could not add to the chain: atoms are lost
}

fn field(line: &str, a: usize, b: usize) -> &str {
    let n = line.len();
    if a >= n {
        return "";
    }
    line.get(a..b.min(n)).unwrap_or("")
}

fn ch(line: &str, i: usize) -> char {
    line.as_bytes().get(i).map(|&c| c as char).unwrap_or('\0')
}

pub fn parse(text: &str) -> Vec<Model> {
    let lines: Vec<&str> = text.split('\n').collect();
    // header: up to the first ATOM/HETATM/MODEL record
    let start = lines
        .iter()
        .position(|l| {
            let r = field(l, 0, 6);
            r == "ATOM  " || r == "HETATM" || r == "MODEL "
        })
        .unwrap_or(lines.len());

    let mut models: Vec<Model> = vec![];
    let mut model_open = false;
    let mut cur_chain: Option<char> = None;
    let mut cur_res: Option<((String, i64, char), String)> = None;
    let mut cur = Cur::None;
    let mut chain_idx = 0usize;

    for raw in &lines[start..] {
        let line = raw.trim_end_matches(['\n', '\r']);
        let rec = field(line, 0, 6);
        if line.trim().is_empty() {
            continue;
        }
        if rec == "ATOM  " || rec == "HETATM" {
            if !model_open {
                models.push(Model::default());
                model_open = true;
            }
            let fullname = field(line, 12, 16).to_string();
            let toks: Vec<&str> = fullname.split_whitespace().collect();
            let mut name = if toks.len() == 1 { toks[0].to_string() } else { fullname.clone() };
            let altloc = ch(line, 16);
            let resname = field(line, 17, 20).trim().to_string();
            let chainid = ch(line, 21);
            let resseq: i64 = match field(line, 22, 26).split_whitespace().next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            let icode = ch(line, 26);
            let hetero = if rec == "HETATM" {
                if resname == "HOH" || resname == "WAT" {
                    "W".to_string()
                } else {
                    format!("H_{}", resname)
                }
            } else {
                " ".to_string()
            };
            let parse_f = |a, b| field(line, a, b).trim().parse::<f64>();
            let (x, y, z) = match (parse_f(30, 38), parse_f(38, 46), parse_f(46, 54)) {
                (Ok(x), Ok(y), Ok(z)) => (x, y, z),
                _ => continue,
            };
            let coord = f32r([x, y, z]);
            let occ = parse_f(54, 60).unwrap_or(0.0);
            let res_id = (hetero.clone(), resseq, icode);
            let model = models.last_mut().unwrap();

            // init_chain / init_residue
            let new_chain = cur_chain != Some(chainid);
            let new_res = new_chain || cur_res.as_ref() != Some(&(res_id.clone(), resname.clone()));
            if new_chain {
                cur_chain = Some(chainid);
                chain_idx = match model.chains.iter().position(|c| c.id == chainid) {
                    Some(i) => i,
                    None => {
                        model.chains.push(Chain { id: chainid, entries: vec![] });
                        model.chains.len() - 1
                    }
                };
            }
            if new_res {
                cur_res = Some((res_id.clone(), resname.clone()));
                cur = init_residue(&mut model.chains[chain_idx], chain_idx, &res_id, &resname, resseq);
            }

            // init_atom
            let res = match cur {
                Cur::Entry(c, e) => {
                    let entry = &mut model.chains[c].entries[e];
                    let s = entry.selected;
                    &mut entry.variants[s]
                }
                _ => continue,
            };
            if let Some(i) = res.find(&name) {
                let dup_full = res.atoms[i].children[res.atoms[i].selected].3.clone();
                if dup_full != fullname {
                    name = fullname.clone();
                }
            }
            let child = (altloc, coord, occ, fullname.clone());
            if altloc != ' ' {
                match res.find(&name) {
                    Some(i) if res.atoms[i].disordered => disordered_add(&mut res.atoms[i], child),
                    Some(i) => {
                        let dup = res.atoms.remove(i);
                        let mut d = Atom {
                            id: name.clone(),
                            children: vec![],
                            selected: 0,
                            last_occ: f64::NEG_INFINITY,
                            disordered: true,
                        };
                        disordered_add(&mut d, child);
                        disordered_add(&mut d, dup.children[dup.selected].clone());
                        res.atoms.push(d);
                    }
                    None => {
                        let mut d = Atom {
                            id: name.clone(),
                            children: vec![],
                            selected: 0,
                            last_occ: f64::NEG_INFINITY,
                            disordered: true,
                        };
                        disordered_add(&mut d, child);
                        res.atoms.push(d);
                    }
                }
            } else if res.find(&name).is_none() {
                res.atoms.push(Atom {
                    id: name,
                    children: vec![child],
                    selected: 0,
                    last_occ: occ,
                    disordered: false,
                });
            } // else: "Atom defined twice" -> ignored
        } else if rec == "MODEL " {
            models.push(Model::default());
            model_open = true;
            cur_chain = None;
            cur_res = None;
        } else if rec == "END   " || rec == "CONECT" {
            break;
        } else if rec == "ENDMDL" {
            model_open = false;
            cur_chain = None;
            cur_res = None;
        }
    }

    // ClaRNA: rename old-style atom names, e.g. O2* -> O2'
    for m in models.iter_mut() {
        for c in m.chains.iter_mut() {
            for e in c.entries.iter_mut() {
                for r in e.variants.iter_mut() {
                    for a in r.atoms.iter_mut() {
                        if is_star_name(&a.id) {
                            a.id = a.id.replace('*', "'");
                        }
                    }
                }
            }
        }
    }
    models
}

/// ^[A-Z]{1,2}[0-9]?\*$
fn is_star_name(s: &str) -> bool {
    let b = s.as_bytes();
    if b.is_empty() || *b.last().unwrap() != b'*' {
        return false;
    }
    let body = &b[..b.len() - 1];
    let n_alpha = body.iter().take_while(|c| c.is_ascii_uppercase()).count();
    if n_alpha < 1 || n_alpha > 2 {
        return false;
    }
    let rest = &body[n_alpha..];
    rest.is_empty() || (rest.len() == 1 && rest[0].is_ascii_digit())
}

fn disordered_add(a: &mut Atom, child: (char, V3, f64, String)) {
    let occ = child.2;
    // Bio: self[altloc] = atom (same altloc replaces)
    let idx = match a.children.iter().position(|c| c.0 == child.0) {
        Some(i) => {
            a.children[i] = child;
            i
        }
        None => {
            a.children.push(child);
            a.children.len() - 1
        }
    };
    if occ > a.last_occ {
        a.last_occ = occ;
        a.selected = idx;
    }
}

fn init_residue(chain: &mut Chain, chain_idx: usize, res_id: &(String, i64, char), resname: &str, resseq: i64) -> Cur {
    let new_res = || Residue { resname: resname.to_string(), resseq, atoms: vec![] };
    if let Some(e) = chain.entries.iter().position(|e| &e.id == res_id) {
        if res_id.0 != " " {
            // Bio: chain.add() fails ("already defined"); atoms go to a residue that is not in the chain
            return Cur::Orphan;
        }
        let entry = &mut chain.entries[e];
        if entry.variants.len() > 1 {
            // disordered residue
            if let Some(v) = entry.variants.iter().position(|r| r.resname == resname) {
                entry.selected = v;
            } else {
                entry.variants.push(new_res());
                entry.selected = entry.variants.len() - 1;
            }
            return Cur::Entry(chain_idx, e);
        }
        if entry.variants[0].resname == resname {
            return Cur::Entry(chain_idx, e);
        }
        let completely_disordered = entry.variants[0].atoms.iter().all(|a| a.altloc() != ' ');
        if !completely_disordered {
            return Cur::None; // "Blank altlocs in duplicate residue": atoms are dropped
        }
        let mut old = chain.entries.remove(e);
        old.variants.push(new_res());
        old.selected = 1;
        chain.entries.push(old);
        return Cur::Entry(chain_idx, chain.entries.len() - 1);
    }
    chain.entries.push(ResEntry { id: res_id.clone(), variants: vec![new_res()], selected: 0 });
    Cur::Entry(chain_idx, chain.entries.len() - 1)
}

/// StrucFile(fn).get_info_chains(): "A 1 76 B 77 90"
pub fn info_chains(text: &str) -> String {
    let mut lines: Vec<String> = vec![];
    for l in text.trim().split('\n') {
        if l.starts_with("ENDMDL") {
            break;
        }
        if l.starts_with("ATOM") || l.starts_with("HETATM") || l.starts_with("TER") || l.starts_with("END") {
            lines.push(l.trim().to_string());
        }
    }
    let mut chains: Vec<(char, String)> = vec![];
    let mut curri: Option<i64> = None;
    let mut chain_prev: Option<char> = None;
    let mut resi_prev: i64 = 0;
    for l in &lines {
        if !(l.starts_with("ATOM") || l.starts_with("HETATM")) {
            continue;
        }
        let resi: i64 = match field(l, 22, 26).trim().parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if curri != Some(resi) {
            let chain_curr = ch(l, 21);
            if let Some(p) = chain_prev {
                if p != chain_curr {
                    let h = &mut chains.iter_mut().find(|c| c.0 == p).unwrap().1;
                    h.push_str(&format!("-{}", resi_prev));
                }
            }
            if !chains.iter().any(|c| c.0 == chain_curr) {
                chains.push((chain_curr, format!("{}:{}", chain_curr, resi)));
            }
            resi_prev = resi;
            chain_prev = Some(chain_curr);
        }
        curri = Some(resi);
    }
    if let Some(p) = chain_prev {
        let h = &mut chains.iter_mut().find(|c| c.0 == p).unwrap().1;
        h.push_str(&format!("-{}", resi_prev));
    }
    let s: Vec<String> = chains.iter().map(|(_, h)| h.replace(':', " ").replace('-', " ")).collect();
    s.join(" ").trim().to_string()
}
