//! clarna_rs: a Rust port of ClaRNA (rna_clarna_run.py, -bp+stack by default).
//!
//!     clarna_rs a.pdb                    # the same output as rna_clarna_run.py -ipdb a.pdb -bp+stack
//!     clarna_rs a.pdb b.pdb ...          # several files (or --batch): each output preceded by "### <file>"
//!     clarna_rs --dump DIR a.pdb ...     # also write DIR/<name>.outCR/.dist.tsv/.contacts.tsv (regression format)
//!
//! Options: -bps | -stack | -bp+stack, -thresh 0.6, --lib DIR (ClaRNAlib/lib)

mod clarna;
mod geom;
mod pdb;
mod program;
mod residue;

use std::io::Write;

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut opts = "bp+stack".to_string();
    let mut min_score = 0.6;
    let mut lib_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../ClaRNAlib/lib").to_string();
    let mut dump: Option<String> = None;
    let mut files = vec![];
    let mut batch = false;
    while !args.is_empty() {
        let a = args.remove(0);
        match a.as_str() {
            "-bps" | "-stack" | "-bp+stack" => opts = a[1..].to_string(),
            "-thresh" => min_score = args.remove(0).parse().expect("-thresh: a number"),
            "--lib" => lib_dir = args.remove(0),
            "--dump" => dump = Some(args.remove(0)),
            "-ipdb" => {}
            "--batch" => batch = true,
            "-h" | "--help" => {
                println!("usage: clarna_rs [-bps|-stack|-bp+stack] [-thresh 0.6] [--lib DIR] [--dump DIR] file.pdb ...");
                return;
            }
            _ => files.push(a),
        }
    }
    let libs: Vec<String> = match opts.as_str() {
        "bps" => vec!["classifier.bp.json.gz"],
        "stack" => vec!["classifier.stacking.json.gz"],
        _ => vec!["classifier.bp.json.gz", "classifier.stacking.json.gz"],
    }
    .iter()
    .map(|l| format!("{}/{}", lib_dir, l))
    .collect();
    let cl = clarna::ClaRNA::new(&libs);

    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    for f in &files {
        if batch || files.len() > 1 {
            writeln!(out, "### {}", f).unwrap();
        }
        let text = match std::fs::read(f) {
            Ok(b) => String::from_utf8_lossy(&b).to_string(),
            Err(e) => {
                writeln!(out, "ERROR: {}: {}", f, e).unwrap();
                continue;
            }
        };
        let models = pdb::parse(&text);
        let st = clarna::extract(&models);
        let g = cl.classify(&st);
        let outcr = format!("Classifier: Clarna\nchains:  {}\n{}", pdb::info_chains(&text), clarna::display(&g, min_score));
        out.write_all(outcr.as_bytes()).unwrap();
        out.flush().unwrap();
        if let Some(dir) = &dump {
            std::fs::create_dir_all(dir).unwrap();
            let name = std::path::Path::new(f).file_stem().unwrap().to_string_lossy().to_string();
            let mut dist = vec![];
            let mut contacts = vec![];
            for (u, v, e) in g.edges() {
                match e {
                    clarna::Edge::Dist { n_type } => dist.push(format!("{}\t{}\t{}\t{:.6}", u, v, n_type, 0.0)),
                    clarna::Edge::Contact(c) => contacts.push(format!(
                        "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}",
                        u,
                        v,
                        c.n_type,
                        c.desc,
                        c.full_desc,
                        c.weight,
                        c.reference.as_deref().unwrap_or("None"),
                        c.is_alt as i32
                    )),
                }
            }
            dist.sort();
            contacts.sort();
            std::fs::write(format!("{}/{}.outCR", dir, name), &outcr).unwrap();
            std::fs::write(format!("{}/{}.dist.tsv", dir, name), dist.join("\n") + "\n").unwrap();
            std::fs::write(format!("{}/{}.contacts.tsv", dir, name), contacts.join("\n") + "\n").unwrap();
        }
    }
}
