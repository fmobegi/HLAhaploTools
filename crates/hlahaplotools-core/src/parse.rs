//! Readers for Hapl-o-Mat's MAC and GLSC input formats.
//!
//! **MAC**: tab-separated, with a header naming each locus twice:
//! ```text
//! id  A      A       B       B  ...
//! 1   26:04  24:DXTW 54:16   07:218
//! ```
//! Typings have no locus prefix. An empty cell (or `-`, `.`, `NA`, `?`) means
//! that allele was not typed.
//!
//! **GLSC**: whitespace-separated; an ID followed by one GL string per locus,
//! in any order:
//! ```text
//! 1   A*11:98+A*30:76   B*51:178+B*44:176/B*07:84   DRB1*13:212+DRB1*08:03P|DRB1*15:02:14+DRB1*04:04:09
//! ```
//! `|` separates alternative genotypes, `+` the two alleles, `/` alternative
//! alleles. A locus left out of a line is treated as missing.
//!
//! Unlike Hapl-o-Mat, missing or half-typed loci do not cause a sample to be
//! dropped: they are kept as [`LocusObs::Missing`] / [`LocusObs::OneAllele`].

use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::Error;
use crate::dataset::{Dataset, Issue, LocusObs, Sample};
use crate::intern::{AlleleId, AlleleTable};
use crate::mac::MacDictionary;
use crate::typing::{TypingKind, classify};

/// Options shared by both readers.
#[derive(Debug, Clone, Default)]
pub struct ParseOptions<'a> {
    /// Loci to keep, in this order. `None` keeps all loci: header order for
    /// MAC files, alphabetical for GLSC files.
    pub loci: Option<Vec<String>>,
    /// Needed to expand MAC codes such as `24:DXTW`. Without it, samples with
    /// MAC codes are skipped and reported as issues.
    pub mac: Option<&'a MacDictionary>,
}

/// Cell values that mean "not typed".
fn is_missing(cell: &str) -> bool {
    matches!(cell.trim(), "" | "-" | "." | "NA" | "N/A" | "?")
}

/// Turns typing strings into allele IDs, expanding MAC codes on the way.
struct Resolver<'a> {
    mac: Option<&'a MacDictionary>,
    table: AlleleTable,
}

impl Resolver<'_> {
    /// Resolve one typing (without locus prefix) to its candidate alleles.
    fn resolve(
        &mut self,
        locus_idx: u16,
        locus: &str,
        typing: &str,
    ) -> Result<Vec<AlleleId>, String> {
        let typing = typing.trim();
        // `NNNN` marks an absent gene (DRB3/4/5 on haplotypes that lack it).
        if typing.eq_ignore_ascii_case("NNNN") {
            return Ok(vec![self.table.intern(locus_idx, &format!("{locus}*NNNN"))]);
        }
        match classify(typing) {
            None => Err(format!("unrecognised typing {locus}*{typing}")),
            Some(TypingKind::Mac) => {
                let (family, code) = typing.split_once(':').expect("MAC typings contain ':'");
                // `01:XX` means "any allele in family 01". NMDP's table has no
                // XX entry (Hapl-o-Mat drops such samples); we read it as the
                // first-field typing `A*01`, which the resolver expands.
                if code == "XX" && self.mac.is_none_or(|d| d.get("XX").is_none()) {
                    return Ok(vec![
                        self.table.intern(locus_idx, &format!("{locus}*{family}")),
                    ]);
                }
                let dict = self.mac.ok_or_else(|| {
                    format!("{locus}*{typing} is a MAC code but no MAC dictionary was given")
                })?;
                let names = dict.expand(locus, typing)?;
                Ok(names
                    .iter()
                    .map(|n| self.table.intern(locus_idx, n))
                    .collect())
            }
            Some(_) => Ok(vec![
                self.table.intern(locus_idx, &format!("{locus}*{typing}")),
            ]),
        }
    }

    /// Resolve a `/`-separated list of full names (`A*01:01/A*01:02`).
    ///
    /// Order and duplicates are kept: each listed entry gets an equal share,
    /// so `X/Y/X` weights X twice, as Hapl-o-Mat does.
    fn resolve_list(
        &mut self,
        locus_idx: u16,
        locus: &str,
        list: &str,
    ) -> Result<Vec<AlleleId>, String> {
        let mut out = Vec::new();
        for full in list.split('/') {
            let typing = strip_locus(locus, full)?;
            out.extend(self.resolve(locus_idx, locus, typing)?);
        }
        Ok(out)
    }
}

/// `A*01:01` -> `01:01`, checking the locus matches.
fn strip_locus<'s>(locus: &str, full: &'s str) -> Result<&'s str, String> {
    match full.split_once('*') {
        Some((l, typing)) if l == locus => Ok(typing),
        Some((l, _)) => Err(format!("allele {full} is not at locus {locus} (found {l})")),
        None => Err(format!("allele {full} has no locus prefix")),
    }
}

fn split_cells(line: &str) -> Vec<&str> {
    if line.contains('\t') {
        line.split('\t').collect()
    } else {
        line.split_whitespace().collect()
    }
}

fn locus_index(i: usize) -> u16 {
    u16::try_from(i).expect("more than 65535 loci")
}

// ------------------------------------------------------- MAC code pre-scan

/// All MAC codes (e.g. `DXTW`) used anywhere in a MAC or GLSC file.
///
/// Pass the result to [`MacDictionary::from_path_only`] to load just those
/// codes instead of the whole ~650 MB NMDP table.
pub fn collect_mac_codes<R: BufRead>(
    reader: R,
) -> Result<std::collections::HashSet<String>, Error> {
    let mut codes = std::collections::HashSet::new();
    for line in reader.lines() {
        let line = line?;
        // Also split on CSV punctuation so wide tables can be scanned too.
        let sep = |c: char| c.is_whitespace() || matches!(c, '+' | '/' | '|' | ',' | ';' | '"');
        for token in line.split(sep) {
            let typing = token.split_once('*').map_or(token, |(_, t)| t);
            if classify(typing) == Some(TypingKind::Mac)
                && let Some((_, code)) = typing.split_once(':')
            {
                codes.insert(code.to_owned());
            }
        }
    }
    Ok(codes)
}

pub fn collect_mac_codes_in_file(
    path: impl AsRef<Path>,
) -> Result<std::collections::HashSet<String>, Error> {
    collect_mac_codes(BufReader::with_capacity(1 << 20, File::open(path)?))
}

// ---------------------------------------------------------------- MAC format

pub fn read_mac<R: BufRead>(reader: R, opts: &ParseOptions) -> Result<Dataset, Error> {
    let mut lines = reader.lines().enumerate();

    // Header: first non-blank line.
    let header = loop {
        match lines.next() {
            Some((_, line)) => {
                let line = line?;
                if !line.trim().is_empty() {
                    break line;
                }
            }
            None => return Err(Error::Format("MAC file is empty".into())),
        }
    };
    let header = header.trim_start_matches('\u{feff}').trim_end_matches('\r');
    let cols = split_cells(header);
    let typing_cols = &cols[1..];
    if typing_cols.is_empty() || !typing_cols.len().is_multiple_of(2) {
        return Err(Error::Format(format!(
            "MAC header must be an ID column followed by two columns per locus, got {} typing columns",
            typing_cols.len()
        )));
    }
    let mut file_loci: Vec<String> = Vec::new();
    for pair in typing_cols.chunks(2) {
        let (a, b) = (pair[0].trim(), pair[1].trim());
        if a != b {
            return Err(Error::Format(format!(
                "MAC header columns {a} and {b} should name the same locus"
            )));
        }
        if file_loci.iter().any(|l| l == a) {
            return Err(Error::Format(format!(
                "locus {a} appears twice in the MAC header"
            )));
        }
        file_loci.push(a.to_owned());
    }

    // Which file columns feed which output locus.
    let loci: Vec<String> = opts.loci.clone().unwrap_or_else(|| file_loci.clone());
    let mut file_pos = Vec::with_capacity(loci.len());
    for locus in &loci {
        let pos = file_loci.iter().position(|l| l == locus).ok_or_else(|| {
            Error::Format(format!(
                "locus {locus} was requested but is not in the MAC header"
            ))
        })?;
        file_pos.push(pos);
    }

    let mut resolver = Resolver {
        mac: opts.mac,
        table: AlleleTable::new(),
    };
    let mut samples = Vec::new();
    let mut issues = Vec::new();

    for (i, line) in lines {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let line_no = i + 1;
        let cells = split_cells(line);
        let id = cells[0].trim().to_owned();
        if cells.len() - 1 != typing_cols.len() {
            issues.push(Issue {
                line: line_no,
                sample_id: Some(id),
                message: format!(
                    "expected {} typings, found {} (missing typings must still be tab-delimited)",
                    typing_cols.len(),
                    cells.len() - 1
                ),
            });
            continue;
        }

        let mut obs = Vec::with_capacity(loci.len());
        let mut failed = None;
        for (k, (locus, &pos)) in loci.iter().zip(&file_pos).enumerate() {
            let (t1, t2) = (cells[1 + 2 * pos], cells[2 + 2 * pos]);
            match mac_locus(&mut resolver, locus_index(k), locus, t1, t2) {
                Ok(o) => obs.push(o),
                Err(msg) => {
                    failed = Some(msg);
                    break;
                }
            }
        }
        match failed {
            None => samples.push(Sample { id, loci: obs }),
            Some(message) => issues.push(Issue {
                line: line_no,
                sample_id: Some(id),
                message,
            }),
        }
    }

    Ok(Dataset {
        loci,
        samples,
        alleles: resolver.table,
        issues,
        notes: Vec::new(),
        copied: Vec::new(),
    })
}

fn mac_locus(
    r: &mut Resolver,
    k: u16,
    locus: &str,
    t1: &str,
    t2: &str,
) -> Result<LocusObs, String> {
    // Tolerate typings written with the locus prefix (`A*01:01`).
    let bare = |t: &'_ str| -> Result<String, String> {
        if t.contains('*') {
            strip_locus(locus, t.trim()).map(str::to_owned)
        } else {
            Ok(t.trim().to_owned())
        }
    };
    Ok(match (is_missing(t1), is_missing(t2)) {
        (true, true) => LocusObs::Missing,
        (false, true) => LocusObs::OneAllele(r.resolve(k, locus, &bare(t1)?)?),
        (true, false) => LocusObs::OneAllele(r.resolve(k, locus, &bare(t2)?)?),
        (false, false) => {
            let a = r.resolve(k, locus, &bare(t1)?)?;
            let b = r.resolve(k, locus, &bare(t2)?)?;
            LocusObs::Genotypes(vec![(a, b)])
        }
    })
}

pub fn read_mac_file(path: impl AsRef<Path>, opts: &ParseOptions) -> Result<Dataset, Error> {
    read_mac(BufReader::with_capacity(1 << 20, File::open(path)?), opts)
}

// --------------------------------------------------------------- GLSC format

pub fn read_glsc<R: BufRead>(reader: R, opts: &ParseOptions) -> Result<Dataset, Error> {
    let mut records = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_start_matches('\u{feff}').trim();
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split_whitespace();
        let id = fields.next().unwrap_or_default().to_owned();
        records.push(GlRecord {
            line: i + 1,
            id,
            fields: fields.map(str::to_owned).collect(),
        });
    }
    Ok(dataset_from_gl(records, Vec::new(), Vec::new(), opts))
}

/// One sample as GL strings, from any GL-based format (GLSC, GLS, HML, NGSengine).
#[derive(Debug, Clone)]
pub(crate) struct GlRecord {
    /// 1-based line in the input (0 when the format has no useful line number).
    pub line: usize,
    pub id: String,
    /// GL strings, one or more loci each (`^` joins loci; `HLA-` prefixes allowed).
    pub fields: Vec<String>,
}

/// Tidy a GL string as written in HML or vendor exports: drop whitespace and
/// `HLA-` prefixes, and split `^`-joined loci into one string per locus.
pub(crate) fn normalise_gl(gl: &str) -> Result<Vec<String>, String> {
    let compact: String = gl.chars().filter(|c| !c.is_whitespace()).collect();
    let compact = compact.replace("HLA-", "");
    if compact.contains('~') {
        return Err(format!(
            "phased GL string {compact} ('~') is not supported; give each locus separately"
        ));
    }
    Ok(compact
        .split('^')
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .collect())
}

/// Build a dataset from GL-string records. `issues` are records already
/// rejected by the reader; `notes` are loci a reader left untyped on purpose.
pub(crate) fn dataset_from_gl(
    records: Vec<GlRecord>,
    mut issues: Vec<Issue>,
    notes: Vec<Issue>,
    opts: &ParseOptions,
) -> Dataset {
    // Normalise every field first, so locus names are known.
    let mut rows: Vec<(usize, String, Result<Vec<String>, String>)> =
        Vec::with_capacity(records.len());
    for r in records {
        let mut fields = Vec::with_capacity(r.fields.len());
        let mut err = None;
        for f in &r.fields {
            match normalise_gl(f) {
                Ok(parts) => fields.extend(parts),
                Err(e) => {
                    err = Some(e);
                    break;
                }
            }
        }
        rows.push((r.line, r.id, err.map_or(Ok(fields), Err)));
    }

    let loci: Vec<String> = match &opts.loci {
        Some(l) => l.clone(),
        None => {
            let mut seen = BTreeSet::new();
            for (_, _, fields) in &rows {
                for field in fields.iter().flatten() {
                    if let Some((locus, _)) = field.split_once('*') {
                        seen.insert(locus.to_owned());
                    }
                }
            }
            seen.into_iter().collect()
        }
    };
    let locus_pos: HashMap<&str, usize> = loci
        .iter()
        .enumerate()
        .map(|(i, l)| (l.as_str(), i))
        .collect();

    let mut resolver = Resolver {
        mac: opts.mac,
        table: AlleleTable::new(),
    };
    let mut samples = Vec::new();

    for (line_no, id, fields) in rows {
        let mut obs = vec![LocusObs::Missing; loci.len()];
        let mut seen = vec![false; loci.len()];

        let result: Result<(), String> = fields.and_then(|fields| {
            fields.iter().try_for_each(|field| {
                let (locus, _) = field
                    .split_once('*')
                    .ok_or_else(|| format!("GL string {field} has no locus prefix"))?;
                let Some(&k) = locus_pos.get(locus) else {
                    return Ok(()); // locus not requested
                };
                if std::mem::replace(&mut seen[k], true) {
                    return Err(format!("locus {locus} appears twice"));
                }
                obs[k] = glsc_locus(&mut resolver, locus_index(k), locus, field)?;
                Ok(())
            })
        });

        match result {
            Ok(()) => samples.push(Sample { id, loci: obs }),
            Err(message) => issues.push(Issue {
                line: line_no,
                sample_id: Some(id),
                message,
            }),
        }
    }

    Dataset {
        loci,
        samples,
        alleles: resolver.table,
        issues,
        notes,
        copied: Vec::new(),
    }
}

// ---------------------------------------------------------------- GLS format

/// Read Hapl-o-Mat's GLS input: a pull file (`sample;glid:glid:...`, one GL ID
/// per locus) and a GLID file (`glid;GL string`). GL ID `0` means the locus
/// was not typed; it is kept as missing, not dropped.
///
/// Hapl-o-Mat's `LOCIORDER` is not needed: each GL string names its locus.
pub fn read_gls<P: BufRead, G: BufRead>(
    pull: P,
    glid: G,
    opts: &ParseOptions,
) -> Result<Dataset, Error> {
    let mut gls: HashMap<u64, String> = HashMap::new();
    for (i, line) in glid.lines().enumerate() {
        let line = line?;
        let line = line.trim_start_matches('\u{feff}').trim();
        if line.is_empty() {
            continue;
        }
        let (n, gl) = line.split_once(';').ok_or_else(|| {
            Error::Format(format!("GLID file line {}: expected `id;GL string`", i + 1))
        })?;
        let n: u64 = n.trim().parse().map_err(|_| {
            Error::Format(format!("GLID file line {}: {n:?} is not a GL ID", i + 1))
        })?;
        gls.insert(n, gl.trim().to_owned());
    }

    let mut records = Vec::new();
    let mut issues = Vec::new();
    for (i, line) in pull.lines().enumerate() {
        let line = line?;
        let line = line.trim_start_matches('\u{feff}').trim();
        if line.is_empty() {
            continue;
        }
        let (id, ids) = line.split_once(';').unwrap_or((line, ""));
        let id = id.trim().to_owned();
        let mut fields = Vec::new();
        let mut err = None;
        for n in ids.split(':').map(str::trim).filter(|s| !s.is_empty()) {
            match n.parse::<u64>() {
                Ok(0) => {} // not typed
                Ok(n) => match gls.get(&n) {
                    Some(gl) => fields.push(gl.clone()),
                    None => {
                        err = Some(format!("GL ID {n} is not in the GLID file"));
                        break;
                    }
                },
                Err(_) => {
                    err = Some(format!("{n:?} is not a GL ID"));
                    break;
                }
            }
        }
        match err {
            None => records.push(GlRecord {
                line: i + 1,
                id,
                fields,
            }),
            Some(message) => issues.push(Issue {
                line: i + 1,
                sample_id: Some(id),
                message,
            }),
        }
    }
    Ok(dataset_from_gl(records, issues, Vec::new(), opts))
}

pub fn read_gls_files(
    pull: impl AsRef<Path>,
    glid: impl AsRef<Path>,
    opts: &ParseOptions,
) -> Result<Dataset, Error> {
    read_gls(
        BufReader::with_capacity(1 << 20, File::open(pull)?),
        BufReader::with_capacity(1 << 20, File::open(glid)?),
        opts,
    )
}

/// Parse one locus GL string, e.g. `A*01:01+A*02:01/A*02:02|A*03:01+A*11:01`.
fn glsc_locus(r: &mut Resolver, k: u16, locus: &str, gl: &str) -> Result<LocusObs, String> {
    // No `+` at all (`DRB4*01:01:01G`, `DRB3*01:01|DRB3*01:02`): one allele
    // reported, as labs often do for DRB3/4/5 carried on one chromosome. The
    // alternatives become candidates for that allele; filling decides the rest.
    if !gl.contains('+') {
        let list = gl.replace('|', "/");
        return Ok(LocusObs::OneAllele(r.resolve_list(k, locus, &list)?));
    }
    let genotypes: Vec<&str> = gl.split('|').collect();
    let mut out = Vec::with_capacity(genotypes.len());
    for genotype in &genotypes {
        let sides: Vec<&str> = genotype.split('+').collect();
        let [left, right] = sides[..] else {
            return Err(format!("genotype {genotype} should have exactly one '+'"));
        };
        match (is_missing(left), is_missing(right)) {
            (false, false) => {
                let a = r.resolve_list(k, locus, left)?;
                let b = r.resolve_list(k, locus, right)?;
                out.push((a, b));
            }
            (true, true) => return Err(format!("genotype {genotype} has no alleles")),
            // Half-typed (`A*01:01+`): only allowed as the sole genotype.
            _ if genotypes.len() == 1 => {
                let known = if is_missing(left) { right } else { left };
                return Ok(LocusObs::OneAllele(r.resolve_list(k, locus, known)?));
            }
            _ => {
                return Err(format!(
                    "half-typed genotype {genotype} inside a genotype list"
                ));
            }
        }
    }
    Ok(LocusObs::Genotypes(out))
}

// ------------------------------------------------------------ wide tables

/// Gene names recognised in column headers like `A1` (no separator). With a
/// separator (`A_1`, `A.1`, `A-1`, `A_ALLELE_1`) these and any loci passed in
/// [`ParseOptions::loci`] are accepted.
const KNOWN_GENES: &[&str] = &[
    "A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V", "DRA", "DRB1", "DRB3", "DRB4", "DRB5",
    "DQA1", "DQB1", "DPA1", "DPB1", "DMA", "DMB", "DOA", "DOB", "MICA", "MICB", "HFE",
];

/// Identifier columns found in a wide table, aligned with [`Dataset::samples`].
#[derive(Debug, Clone, Default)]
pub struct TableIds {
    /// Value of the family column (`FAMILY_ID`, `Family`, ...) per sample, if present.
    pub family: Option<Vec<String>>,
    /// Value of the member column (`Family_Member`, `Relationship`, ...) per sample.
    pub member: Option<Vec<String>>,
}

/// `A_1` -> ("A", 0), `DRB1.2` -> ("DRB1", 1), `HLA-A_ALLELE_1` -> ("A", 0).
fn allele_column(header: &str, extra_genes: &[String]) -> Option<(String, usize)> {
    let mut h = header.trim().to_ascii_uppercase();
    for prefix in ["HLA-", "HLA_", "HLA."] {
        if let Some(rest) = h.strip_prefix(prefix) {
            h = rest.to_string();
        }
    }
    let which = match h.chars().last()? {
        '1' => 0,
        '2' => 1,
        _ => return None,
    };
    let mut rest = &h[..h.len() - 1];
    let is_sep = |c: char| matches!(c, '_' | '.' | '-' | ' ');
    rest = rest.strip_suffix(is_sep).unwrap_or(rest);
    if let Some(r) = rest.strip_suffix("ALLELE") {
        rest = r.trim_end_matches(is_sep);
    }
    // Only known genes (or loci the caller asked for) count, so that columns
    // like "Batch_1" are not mistaken for alleles.
    let known =
        KNOWN_GENES.contains(&rest) || extra_genes.iter().any(|g| g.eq_ignore_ascii_case(rest));
    known.then(|| (rest.to_string(), which))
}

/// Split one line on `delim`, honouring double quotes.
fn split_quoted(line: &str, delim: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut field = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' if in_quotes && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => in_quotes = !in_quotes,
            c if c == delim && !in_quotes => out.push(std::mem::take(&mut field)),
            c => field.push(c),
        }
    }
    out.push(field);
    out.into_iter().map(|f| f.trim().to_string()).collect()
}

/// Cell values in wide tables that mean "not typed".
fn is_missing_cell(cell: &str) -> bool {
    let c = cell.trim();
    c.is_empty()
        || ["na", "n/a", "null", "none", "unknown", "-", ".", "?", "nt"]
            .iter()
            .any(|m| c.eq_ignore_ascii_case(m))
}

/// Normalise one cell to a `/`-list of full names: `01:01`, `*01:01`,
/// `A*01:01` and `HLA-A*01:01` all become `A*01:01`.
fn normalise_cell(gene: &str, cell: &str) -> Result<String, String> {
    let mut parts = Vec::new();
    for part in cell.split('/') {
        let p = part.trim();
        let p = p.strip_prefix("HLA-").unwrap_or(p);
        let typing = match p.split_once('*') {
            Some(("", t)) => t,
            Some((l, t)) if l.eq_ignore_ascii_case(gene) => t,
            Some((l, _)) => return Err(format!("allele {p} found in a {gene} column (locus {l})")),
            None => p,
        };
        if typing.is_empty() {
            return Err(format!("empty allele in {gene} cell {cell:?}"));
        }
        parts.push(format!("{gene}*{typing}"));
    }
    Ok(parts.join("/"))
}

/// Read a wide table: one row per person, two columns per locus.
///
/// ```text
/// FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2,...
/// fam01,F,A*01:01,A*02:01,B*08:01,B*07:02,...
/// ```
///
/// The delimiter (tab, comma or semicolon) is taken from the header. Family
/// and member columns are returned in [`TableIds`]; a `Sample_ID`-like column
/// gives the sample IDs, otherwise `family/member` or the row number is used.
pub fn read_table<R: BufRead>(
    reader: R,
    opts: &ParseOptions,
) -> Result<(Dataset, TableIds), Error> {
    let mut lines = reader.lines().enumerate();
    let header = loop {
        match lines.next() {
            Some((_, line)) => {
                let line = line?;
                if !line.trim().is_empty() {
                    break line;
                }
            }
            None => return Err(Error::Format("table is empty".into())),
        }
    };
    let header = header.trim_start_matches('\u{feff}').trim_end_matches('\r');
    let delim = if header.contains('\t') {
        '\t'
    } else if header.contains(',') {
        ','
    } else {
        ';'
    };
    let cols = split_quoted(header, delim);
    let find = |names: &[&str]| {
        cols.iter().position(|c| {
            let c = c.to_ascii_lowercase().replace([' ', '-'], "_");
            names.contains(&c.as_str())
        })
    };
    let family_col = find(&[
        "family_id",
        "familyid",
        "family",
        "fam_id",
        "famid",
        "family_no",
    ]);
    let member_col = find(&[
        "family_member",
        "familymember",
        "member",
        "relationship",
        "relation",
        "role",
    ]);
    let id_col = find(&[
        "sample_id",
        "sampleid",
        "sample",
        "id",
        "patient_id",
        "subject_id",
        "donor_id",
    ]);

    let extra = opts.loci.clone().unwrap_or_default();
    let mut file_loci: Vec<String> = Vec::new();
    let mut positions: Vec<[Option<usize>; 2]> = Vec::new();
    for (i, c) in cols.iter().enumerate() {
        if Some(i) == family_col || Some(i) == member_col || Some(i) == id_col {
            continue;
        }
        if let Some((gene, which)) = allele_column(c, &extra) {
            let k = match file_loci.iter().position(|g| *g == gene) {
                Some(k) => k,
                None => {
                    file_loci.push(gene);
                    positions.push([None, None]);
                    file_loci.len() - 1
                }
            };
            if positions[k][which].replace(i).is_some() {
                return Err(Error::Format(format!(
                    "two columns for {}_{}",
                    file_loci[k],
                    which + 1
                )));
            }
        }
    }
    for (gene, pos) in file_loci.iter().zip(&positions) {
        if pos[0].is_none() || pos[1].is_none() {
            return Err(Error::Format(format!(
                "locus {gene} needs two columns ({gene}_1 and {gene}_2)"
            )));
        }
    }
    if file_loci.is_empty() {
        return Err(Error::Format(
            "no allele columns found (expected names like A_1, A_2, DRB1_1, ...)".into(),
        ));
    }
    let loci: Vec<String> = opts.loci.clone().unwrap_or_else(|| file_loci.clone());
    let mut use_cols = Vec::with_capacity(loci.len());
    for locus in &loci {
        let k = file_loci
            .iter()
            .position(|g| g.eq_ignore_ascii_case(locus))
            .ok_or_else(|| {
                Error::Format(format!("locus {locus} was requested but has no columns"))
            })?;
        use_cols.push([positions[k][0].unwrap(), positions[k][1].unwrap()]);
    }

    let mut resolver = Resolver {
        mac: opts.mac,
        table: AlleleTable::new(),
    };
    let (mut samples, mut issues) = (Vec::new(), Vec::new());
    let (mut fams, mut members) = (Vec::new(), Vec::new());
    for (i, line) in lines {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let cells = split_quoted(line, delim);
        let cell = |i: Option<usize>| {
            i.and_then(|i| cells.get(i))
                .map(String::as_str)
                .unwrap_or("")
        };
        let (fam, member) = (cell(family_col).to_string(), cell(member_col).to_string());
        let id = match id_col {
            Some(_) if !cell(id_col).is_empty() => cell(id_col).to_string(),
            _ if family_col.is_some() && member_col.is_some() => format!("{fam}/{member}"),
            _ => format!("row{}", i + 1),
        };
        let mut obs = Vec::with_capacity(loci.len());
        let mut failed = None;
        for (k, (locus, [c1, c2])) in loci.iter().zip(&use_cols).enumerate() {
            let (t1, t2) = (cell(Some(*c1)), cell(Some(*c2)));
            let side = |r: &mut Resolver, t: &str| -> Result<Vec<AlleleId>, String> {
                r.resolve_list(locus_index(k), locus, &normalise_cell(locus, t)?)
            };
            let result = match (is_missing_cell(t1), is_missing_cell(t2)) {
                (true, true) => Ok(LocusObs::Missing),
                (false, true) => side(&mut resolver, t1).map(LocusObs::OneAllele),
                (true, false) => side(&mut resolver, t2).map(LocusObs::OneAllele),
                (false, false) => side(&mut resolver, t1)
                    .and_then(|a| Ok(LocusObs::Genotypes(vec![(a, side(&mut resolver, t2)?)]))),
            };
            match result {
                Ok(o) => obs.push(o),
                Err(e) => {
                    failed = Some(e);
                    break;
                }
            }
        }
        match failed {
            None => {
                samples.push(Sample { id, loci: obs });
                fams.push(fam);
                members.push(member);
            }
            Some(message) => issues.push(Issue {
                line: i + 1,
                sample_id: Some(id),
                message,
            }),
        }
    }
    let ids = TableIds {
        family: family_col.map(|_| fams),
        member: member_col.map(|_| members),
    };
    Ok((
        Dataset {
            loci,
            samples,
            alleles: resolver.table,
            issues,
            notes: Vec::new(),
            copied: Vec::new(),
        },
        ids,
    ))
}

pub fn read_table_file(
    path: impl AsRef<Path>,
    opts: &ParseOptions,
) -> Result<(Dataset, TableIds), Error> {
    read_table(BufReader::with_capacity(1 << 20, File::open(path)?), opts)
}

pub fn read_glsc_file(path: impl AsRef<Path>, opts: &ParseOptions) -> Result<Dataset, Error> {
    read_glsc(BufReader::with_capacity(1 << 20, File::open(path)?), opts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn macs() -> MacDictionary {
        let mut d = MacDictionary::new();
        d.insert("AB", "01/02");
        d
    }

    fn gl(ds: &Dataset, sample: usize, locus: usize) -> String {
        ds.samples[sample].loci[locus].to_gl(&ds.alleles)
    }

    #[test]
    fn mac_file_with_codes_and_missing_typings() {
        let input = "id\tA\tA\tB\tB\n\
                     s1\t01:01\t02:AB\t07:02\t08:01\n\
                     s2\t\t03:01\t\t\n\
                     s3\t01:01\t99:ZZ\t07:02\t08:01\n\
                     s4\t01:XX\t02:01\t07:02\t08:01\n";
        let dict = macs();
        let ds = read_mac(
            input.as_bytes(),
            &ParseOptions {
                loci: None,
                mac: Some(&dict),
            },
        )
        .unwrap();

        assert_eq!(ds.loci, ["A", "B"]);
        assert_eq!(ds.samples.len(), 3);
        assert_eq!(gl(&ds, 0, 0), "A*01:01+A*02:01/A*02:02");
        assert_eq!(gl(&ds, 0, 1), "B*07:02+B*08:01");
        assert_eq!(gl(&ds, 1, 0), "A*03:01+");
        assert!(ds.samples[1].loci[1].is_missing());
        // 01:XX is read as the first-field typing A*01.
        assert_eq!(gl(&ds, 2, 0), "A*01+A*02:01");

        // Unknown MAC code -> sample skipped and reported.
        assert_eq!(ds.issues.len(), 1);
        assert_eq!(ds.issues[0].line, 4);
        assert!(ds.issues[0].message.contains("ZZ"));
        assert_eq!(ds.missing_counts(), [0, 1]);
        assert_eq!(ds.one_allele_counts(), [1, 0]);
    }

    #[test]
    fn mac_locus_selection_and_bad_header() {
        let input = "id\tA\tA\tB\tB\ns1\t01:01\t02:01\t07:02\t08:01\n";
        let opts = ParseOptions {
            loci: Some(vec!["B".into()]),
            mac: None,
        };
        let ds = read_mac(input.as_bytes(), &opts).unwrap();
        assert_eq!(ds.loci, ["B"]);
        assert_eq!(gl(&ds, 0, 0), "B*07:02+B*08:01");

        let opts = ParseOptions {
            loci: Some(vec!["C".into()]),
            mac: None,
        };
        assert!(read_mac(input.as_bytes(), &opts).is_err());
        assert!(read_mac("id\tA\tB\n".as_bytes(), &ParseOptions::default()).is_err());
    }

    #[test]
    fn glsc_lists_ambiguities_and_missing_loci() {
        let input = "1\tB*51:01+B*44:02/B*07:02\tA*01:01+A*02:01\n\
                     2\tA*03:01+A*11:01|A*03:02+A*11:01\n\
                     3\tA*01:01+A*01:AB\tB*08:01+\n\
                     4\tA*01:01+B*08:01\n";
        let dict = macs();
        let ds = read_glsc(
            input.as_bytes(),
            &ParseOptions {
                loci: None,
                mac: Some(&dict),
            },
        )
        .unwrap();

        assert_eq!(ds.loci, ["A", "B"]);
        assert_eq!(ds.samples.len(), 3);
        assert_eq!(gl(&ds, 0, 0), "A*01:01+A*02:01");
        assert_eq!(gl(&ds, 0, 1), "B*51:01+B*44:02/B*07:02");
        assert_eq!(gl(&ds, 1, 0), "A*03:01+A*11:01|A*03:02+A*11:01");
        assert!(ds.samples[1].loci[1].is_missing());
        // A*01:AB expands to 01:01/01:02.
        assert_eq!(gl(&ds, 2, 0), "A*01:01+A*01:01/A*01:02");
        assert_eq!(gl(&ds, 2, 1), "B*08:01+");

        assert_eq!(ds.issues.len(), 1);
        assert!(ds.issues[0].message.contains("not at locus A"));
    }

    #[test]
    fn wide_table_with_family_columns_and_mixed_notation() {
        let input = "Family ID,Member,HLA-A_1,A.2,B1,B_ALLELE_2,Batch_1\n\
                     fam1,F,A*01:01,01:01,*08:01,B*07:02,x\n\
                     fam1,M,02:01/02:02,NA,08:01,08:01,y\n\
                     fam1,C1,\"A*01:01\",02:01,,,z\n\
                     fam1,C2,B*01:01,02:01,08:01,08:01,z\n";
        let (ds, ids) = read_table(input.as_bytes(), &ParseOptions::default()).unwrap();
        assert_eq!(ds.loci, ["A", "B"]); // Batch_1 is not an allele column
        assert_eq!(ds.samples.len(), 3); // C2 has a B allele in an A column
        assert_eq!(ds.issues.len(), 1);
        assert!(ds.issues[0].message.contains("B*01:01"));
        assert_eq!(ds.samples[0].id, "fam1/F");
        assert_eq!(gl(&ds, 0, 0), "A*01:01+A*01:01");
        assert_eq!(gl(&ds, 0, 1), "B*08:01+B*07:02");
        assert_eq!(gl(&ds, 1, 0), "A*02:01/A*02:02+");
        assert!(ds.samples[2].loci[1].is_missing());
        assert_eq!(ids.family.unwrap(), ["fam1", "fam1", "fam1"]);
        assert_eq!(ids.member.unwrap(), ["F", "M", "C1"]);
    }

    #[test]
    fn wide_table_tab_separated_with_sample_ids() {
        let input = "SampleID\tA_1\tA_2\ns1\t01:01\t02:01\n";
        let (ds, ids) = read_table(input.as_bytes(), &ParseOptions::default()).unwrap();
        assert_eq!(ds.samples[0].id, "s1");
        assert!(ids.family.is_none() && ids.member.is_none());
        assert!(read_table("id\tA_1\n".as_bytes(), &ParseOptions::default()).is_err());
    }

    #[test]
    fn collects_mac_codes_from_both_formats() {
        let mac = "id\tA\tA\tDPB1\tDPB1\n1\t24:DXTW\t02:01\t09:EMP\t04:01\n";
        let glsc = "1\tA*01:01+A*01:AB/A*02:CD|A*03:01+A*03:01\tB*15:03:01G+B*08:01P\n";
        let got = |s: &str| {
            let mut v: Vec<String> = collect_mac_codes(s.as_bytes())
                .unwrap()
                .into_iter()
                .collect();
            v.sort();
            v
        };
        assert_eq!(got(mac), ["DXTW", "EMP"]);
        assert_eq!(got(glsc), ["AB", "CD"]);
    }

    #[test]
    fn glsc_requested_loci_absent_from_line_are_missing() {
        let input = "1\tA*01:01+A*02:01\n";
        let opts = ParseOptions {
            loci: Some(vec!["A".into(), "DRB1".into()]),
            mac: None,
        };
        let ds = read_glsc(input.as_bytes(), &opts).unwrap();
        assert_eq!(ds.samples[0].loci[1], LocusObs::Missing);
    }

    #[test]
    fn gls_pull_and_glid_with_untyped_loci() {
        let glid = "1;A*01:01+A*02:01\n2;B*08:01+B*07:02/B*07:03\n3;HLA-A*03:01+HLA-A*11:01\n";
        let pull = "s1;1:2\ns2;3:0\ns3;1:9\n";
        let ds = read_gls(pull.as_bytes(), glid.as_bytes(), &ParseOptions::default()).unwrap();
        assert_eq!(ds.loci, ["A", "B"]);
        assert_eq!(ds.samples.len(), 2);
        assert_eq!(gl(&ds, 0, 1), "B*08:01+B*07:02/B*07:03");
        // GL ID 0: untyped, not dropped; HLA- prefixes are accepted.
        assert_eq!(gl(&ds, 1, 0), "A*03:01+A*11:01");
        assert!(ds.samples[1].loci[1].is_missing());
        assert_eq!(ds.issues.len(), 1);
        assert!(ds.issues[0].message.contains("GL ID 9"));
    }

    #[test]
    fn single_allele_without_plus_is_one_allele_typed() {
        let input = "s1 DRB4*01:01:01G DRB3*01:01|DRB3*01:02\n";
        let ds = read_glsc(input.as_bytes(), &ParseOptions::default()).unwrap();
        assert_eq!(ds.loci, ["DRB3", "DRB4"]);
        assert_eq!(gl(&ds, 0, 0), "DRB3*01:01/DRB3*01:02+");
        assert_eq!(gl(&ds, 0, 1), "DRB4*01:01:01G+");
    }

    #[test]
    fn glsc_accepts_hla_prefix_and_joined_loci() {
        let input = "s1 HLA-A*01:01+HLA-A*02:01^HLA-B*08:01+HLA-B*07:02\n";
        let ds = read_glsc(input.as_bytes(), &ParseOptions::default()).unwrap();
        assert_eq!(ds.loci, ["A", "B"]);
        assert_eq!(gl(&ds, 0, 1), "B*08:01+B*07:02");
    }
}
