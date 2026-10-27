//! Family (pedigree) data: reading, Mendelian checks, phasing, and estimation.
//!
//! Relatives share haplotypes, so treating family members as unrelated
//! samples counts the same haplotypes several times. Here:
//!
//! * [`read_family_table`] reads a wide table with family and member columns
//!   (`FAMILY_ID`, `Family_Member` = F / M / C1, C2, ...).
//! * [`check_mendelian`] flags children whose typing cannot come from their
//!   parents (typing errors, sample swaps, non-paternity).
//! * [`phase_families`] works out the parents' four haplotypes from the
//!   children: it tries every way the children could have inherited them
//!   and keeps the ways that fit the most loci. Untyped parents are rebuilt
//!   from the children as far as the data allows.
//! * [`estimate_families`] estimates frequencies from founders (parents)
//!   only, using the phase the family revealed.
//! * [`compare_to_families`] checks an EM estimate against family phasing.
//!
//! Everything works at the target resolution (g, 2f, ...), so a parent typed
//! at 2 fields and a child typed at 3 are compared on equal terms.

use std::collections::BTreeSet;
use std::path::Path;

use rustc_hash::FxHashMap;

use crate::Error;
use crate::dataset::{Dataset, Issue, LocusObs};
use crate::estimate::{
    EmOptions, Estimate, LocusGeno, Resolved, ResolvedSample, WILD, estimate_resolved,
    resolve_dataset,
};
use crate::intern::{AlleleId, AlleleTable};
use crate::nomenclature::Nomenclature;
use crate::parse::{ParseOptions, TableIds, read_table, read_table_file};
use crate::resolve::Resolution;

/// Children used in the phasing search; the search grows as 4^(children - 1).
pub const MAX_PHASING_CHILDREN: usize = 8;

// ------------------------------------------------------------ pedigree

/// A person's role in a nuclear family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    Father,
    Mother,
    Child,
}

/// Read a role label: F / Father / Dad / P, M / Mother / Mum, C1 / Child2 / Sib3 ...
pub fn parse_role(label: &str) -> Option<Role> {
    let l = label.trim().to_ascii_lowercase();
    let base = l
        .trim_end_matches(|c: char| c.is_ascii_digit())
        .trim_end_matches(['_', '-', ' ', '.']);
    match base {
        "f" | "fa" | "father" | "dad" | "p" | "pat" | "paternal" => Some(Role::Father),
        "m" | "mo" | "mother" | "mom" | "mum" | "mat" | "maternal" => Some(Role::Mother),
        "c" | "ch" | "child" | "children" | "kid" | "s" | "sib" | "sibling" | "offspring"
        | "son" | "daughter" | "proband" | "index" => Some(Role::Child),
        _ => None,
    }
}

/// One nuclear family; members are indices into [`Dataset::samples`].
#[derive(Debug, Clone, Default)]
pub struct Family {
    pub id: String,
    pub father: Option<usize>,
    pub mother: Option<usize>,
    pub children: Vec<usize>,
}

/// Families, plus samples not in any family.
#[derive(Debug, Clone, Default)]
pub struct Pedigree {
    pub families: Vec<Family>,
    /// Samples with no family ID: used as unrelated individuals.
    pub unrelated: Vec<usize>,
    /// Members left out (unknown role, duplicate parent, ...).
    pub issues: Vec<Issue>,
}

impl Pedigree {
    /// Build from the family and member columns of a table.
    pub fn from_ids(ds: &Dataset, ids: &TableIds) -> Pedigree {
        let mut ped = Pedigree::default();
        let Some(families) = &ids.family else {
            ped.unrelated = (0..ds.samples.len()).collect();
            return ped;
        };
        let mut index: FxHashMap<&str, usize> = FxHashMap::default();
        for (i, fam) in families.iter().enumerate() {
            let sample_id = ds.samples[i].id.clone();
            let issue = |message: String| Issue {
                line: 0,
                sample_id: Some(sample_id.clone()),
                message,
            };
            if fam.trim().is_empty() {
                ped.unrelated.push(i);
                continue;
            }
            let label = ids.member.as_ref().map_or("", |m| m[i].as_str());
            let Some(role) = parse_role(label) else {
                ped.issues.push(issue(format!(
                    "family {fam}: unknown member role {label:?}; left out (use F, M, C1, C2, ...)"
                )));
                continue;
            };
            let k = *index.entry(fam.as_str()).or_insert_with(|| {
                ped.families.push(Family {
                    id: fam.clone(),
                    ..Default::default()
                });
                ped.families.len() - 1
            });
            let family = &mut ped.families[k];
            let slot = match role {
                Role::Father => &mut family.father,
                Role::Mother => &mut family.mother,
                Role::Child => {
                    family.children.push(i);
                    continue;
                }
            };
            if slot.is_some() {
                ped.issues.push(issue(format!(
                    "family {fam}: second {role:?} ({label}); left out"
                )));
            } else {
                *slot = Some(i);
            }
        }
        // Children in birth-label order (C, C1, C2, ..., C10), whatever the row
        // order: the lowest-numbered child defines haplotypes A and C.
        if let Some(members) = &ids.member {
            for family in &mut ped.families {
                family.children.sort_by_key(|&c| child_number(&members[c]));
            }
        }
        ped
    }
}

/// The number in a child label (`C2` -> 2, `Child10` -> 10); 0 when there is none.
pub fn child_number(label: &str) -> u64 {
    let digits: String = label
        .trim()
        .chars()
        .rev()
        .take_while(char::is_ascii_digit)
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
        .collect();
    digits.parse().unwrap_or(0)
}

/// Labels of a parent's two haplotypes: A/B for the father, C/D for the mother.
pub const HAPLOTYPE_LABELS: [[char; 2]; 2] = [['A', 'B'], ['C', 'D']];

/// Read a family table (see [`crate::parse::read_table`] for the layout).
pub fn read_family_table(
    path: impl AsRef<Path>,
    opts: &ParseOptions,
) -> Result<(Dataset, Pedigree), Error> {
    let (ds, ids) = read_table_file(path, opts)?;
    let ped = Pedigree::from_ids(&ds, &ids);
    Ok((ds, ped))
}

/// Same as [`read_family_table`], from any reader.
pub fn read_family_table_from<R: std::io::BufRead>(
    reader: R,
    opts: &ParseOptions,
) -> Result<(Dataset, Pedigree), Error> {
    let (ds, ids) = read_table(reader, opts)?;
    let ped = Pedigree::from_ids(&ds, &ids);
    Ok((ds, ped))
}

// ------------------------------------------------------------ Mendelian check

/// A problem found in one family.
#[derive(Debug, Clone, PartialEq)]
pub struct MendelIssue {
    pub family: String,
    pub sample_id: Option<String>,
    pub locus: Option<String>,
    pub message: String,
}

pub(crate) fn geno_text(g: &LocusGeno, codes: &AlleleTable) -> String {
    let name = |c: u32| {
        if c == WILD {
            "?".to_string()
        } else {
            codes.name(AlleleId(c)).to_string()
        }
    };
    format!("{}+{}", name(g.a), name(g.b))
}

pub(crate) fn options_text(opts: &[(LocusGeno, f64)], codes: &AlleleTable) -> String {
    opts.iter()
        .map(|(g, _)| geno_text(g, codes))
        .collect::<Vec<_>>()
        .join(" | ")
}

/// Could a parent with genotype `g` pass on allele `x`? Unknown alleles could be anything.
pub(crate) fn parent_has(g: &LocusGeno, x: u32) -> bool {
    g.a == WILD || g.b == WILD || g.a == x || g.b == x
}

/// Is child genotype `c` possible from parents `f` and `m`?
pub(crate) fn child_possible(c: &LocusGeno, f: &LocusGeno, m: &LocusGeno) -> bool {
    match (c.a, c.b) {
        (WILD, WILD) => true,
        (x, WILD) => parent_has(f, x) || parent_has(m, x),
        (x, y) => (parent_has(f, x) && parent_has(m, y)) || (parent_has(f, y) && parent_has(m, x)),
    }
}

pub(crate) const MISSING_OPTS: &[(LocusGeno, f64)] = &[(LocusGeno::MISSING, 1.0)];

fn mendel_resolved(
    ds: &Dataset,
    ped: &Pedigree,
    resolved: &[Option<ResolvedSample>],
    codes: &AlleleTable,
) -> Vec<MendelIssue> {
    let mut out = Vec::new();
    for fam in &ped.families {
        let issue = |sample: Option<usize>, locus: Option<usize>, message: String| MendelIssue {
            family: fam.id.clone(),
            sample_id: sample.map(|s| ds.samples[s].id.clone()),
            locus: locus.map(|l| ds.loci[l].clone()),
            message,
        };
        let typed = |i: Option<usize>| i.and_then(|i| resolved[i].as_ref());
        let (father, mother) = (typed(fam.father), typed(fam.mother));
        let n_children = fam
            .children
            .iter()
            .filter(|&&c| resolved[c].is_some())
            .count();
        if n_children == 0 {
            out.push(issue(
                None,
                None,
                "no typed children: cannot be phased".into(),
            ));
        }
        if father.is_none() && mother.is_none() {
            if n_children < 2 {
                out.push(issue(
                    None,
                    None,
                    "no typed parents and fewer than 2 children: cannot be phased".into(),
                ));
            }
            continue; // nothing to check children against
        }
        for &c in &fam.children {
            let Some(child) = resolved[c].as_ref() else {
                continue;
            };
            for l in 0..ds.loci.len() {
                let f = father.map_or(MISSING_OPTS, |p| p.per_locus[l].as_slice());
                let m = mother.map_or(MISSING_OPTS, |p| p.per_locus[l].as_slice());
                let ok = child.per_locus[l].iter().any(|(cg, _)| {
                    f.iter()
                        .any(|(fg, _)| m.iter().any(|(mg, _)| child_possible(cg, fg, mg)))
                });
                if !ok {
                    out.push(issue(
                        Some(c),
                        Some(l),
                        format!(
                            "child {} not possible from father {} and mother {}",
                            options_text(&child.per_locus[l], codes),
                            options_text(f, codes),
                            options_text(m, codes)
                        ),
                    ));
                }
            }
        }
    }
    out
}

/// Most members whose copied alleles are tried in every combination; above
/// this, all copies at the locus are undone together.
const MAX_COPIES_SEARCHED: usize = 12;

/// Resolve the family data, then undo homozygous copies that break Mendelian
/// inheritance.
///
/// [`crate::fill_dataset`] reads a blank second allele as homozygous. In a
/// family that reading can be checked: when no inheritance pattern fits a
/// locus, the copies made there are undone (the blank is read as untyped) for
/// the fewest members that make the locus fit. When several smallest sets
/// work, all of their members are undone, since none can be preferred. When
/// even undoing every copy does not help, the copies stay and the locus is
/// reported as inconsistent as usual. Returns one list of notes per family.
pub(crate) fn resolve_families<'a>(
    ds: &'a Dataset,
    ped: &Pedigree,
    nom: &'a Nomenclature,
    resolutions: &'a [Resolution],
) -> Result<(Resolved<'a>, Vec<Vec<String>>), Error> {
    let (mut resolver, mut resolved, issues) = resolve_dataset(ds, nom, resolutions)?;
    let mut notes = vec![Vec::new(); ped.families.len()];
    let copied: FxHashMap<(usize, usize), &[AlleleId]> = ds
        .copied
        .iter()
        .map(|c| ((c.sample, c.locus), c.typed.as_slice()))
        .collect();
    if copied.is_empty() {
        return Ok(((resolver, resolved, issues), notes));
    }
    for (fi, fam) in ped.families.iter().enumerate() {
        let children: Vec<usize> = fam
            .children
            .iter()
            .copied()
            .filter(|&c| resolved[c].is_some())
            .collect();
        if children.is_empty() {
            continue;
        }
        let members: Vec<usize> = fam
            .father
            .into_iter()
            .chain(fam.mother)
            .chain(children.iter().copied())
            .filter(|&m| resolved[m].is_some())
            .collect();
        for k in 0..ds.loci.len() {
            // Members copied at this locus, with their genotypes read as half-typed.
            let mut undo: Vec<(usize, Vec<(LocusGeno, f64)>)> = Vec::new();
            for &m in &members {
                if let Some(typed) = copied.get(&(m, k))
                    && let Ok(half) = resolver.locus(&LocusObs::OneAllele(typed.to_vec()))
                {
                    undo.push((m, half));
                }
            }
            if undo.is_empty() {
                continue;
            }
            // Does the locus fit when the members in `mask` are read as half-typed?
            let fits = |mask: u32| -> bool {
                let get = |m: usize| -> &[(LocusGeno, f64)] {
                    match undo.iter().position(|(u, _)| *u == m) {
                        Some(j) if mask & (1 << j) != 0 => &undo[j].1,
                        _ => &resolved[m].as_ref().unwrap().per_locus[k],
                    }
                };
                let parent = |p: Option<usize>| p.filter(|&i| resolved[i].is_some()).map(get);
                let kids: Vec<&[(LocusGeno, f64)]> = children.iter().map(|&c| get(c)).collect();
                LocusModel::build(parent(fam.father), parent(fam.mother), &kids).fits()
            };
            if fits(0) {
                continue;
            }
            let n = undo.len();
            let chosen: u32 = if n > MAX_COPIES_SEARCHED {
                let all = (1u32 << n.min(31)) - 1;
                if fits(all) { all } else { 0 }
            } else {
                (1..=n as u32)
                    .map(|size| {
                        (1u32..1 << n)
                            .filter(|m| m.count_ones() == size && fits(*m))
                            .fold(0, |acc, m| acc | m)
                    })
                    .find(|&union| union != 0)
                    .unwrap_or(0)
            };
            if chosen == 0 {
                continue;
            }
            let mut ids = Vec::new();
            for (j, (m, half)) in undo.into_iter().enumerate() {
                if chosen & (1 << j) != 0 {
                    ids.push(ds.samples[m].id.clone());
                    resolved[m].as_mut().unwrap().per_locus[k] = half;
                }
            }
            notes[fi].push(format!(
                "{}: blank second allele of {} read as untyped, not homozygous (copying it breaks Mendelian inheritance)",
                ds.loci[k],
                ids.join(", ")
            ));
        }
    }
    Ok(((resolver, resolved, issues), notes))
}

/// Check every child against its parents, locus by locus, at the target
/// resolution. Also reports families that cannot be phased.
pub fn check_mendelian(
    ds: &Dataset,
    ped: &Pedigree,
    nom: &Nomenclature,
    resolutions: &[Resolution],
) -> Result<Vec<MendelIssue>, Error> {
    let ((resolver, resolved, _), _) = resolve_families(ds, ped, nom, resolutions)?;
    Ok(mendel_resolved(ds, ped, &resolved, &resolver.codes))
}

// ------------------------------------------------------------ phasing

/// What phasing found for one parent at one locus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ParentLocus {
    /// Allele on haplotype 1 and on haplotype 2 (WILD if not known).
    Phased(u32, u32),
    /// Genotype known from the parent's own typing, but not which haplotype each allele is on.
    Unphased,
    /// Nothing known.
    Unknown,
}

/// Overall outcome of phasing one family.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PhaseStatus {
    /// All four parental haplotypes known at every locus.
    Complete,
    /// Some loci ambiguous, untyped or inconsistent (marked `LOCUS*?`).
    Partial,
    /// Could not be phased; the reason is given.
    Failed(String),
}

/// Which parental haplotype a child received from each parent: `A` or `B`
/// from the father, `C` or `D` from the mother.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChildInheritance {
    pub sample_id: String,
    /// `None` if the data does not decide it.
    pub paternal: Option<char>,
    pub maternal: Option<char>,
}

impl ChildInheritance {
    /// Both labels as text, e.g. `"AC"`, or `"B?"` when the maternal one is undecided.
    pub fn labels(&self) -> String {
        [self.paternal, self.maternal]
            .iter()
            .map(|c| c.unwrap_or('?'))
            .collect()
    }
}

/// Phasing result for one family.
#[derive(Debug, Clone)]
pub struct FamilyPhase {
    pub family: String,
    /// Sample IDs of the typed parents (`None` if not typed).
    pub father: Option<String>,
    pub mother: Option<String>,
    /// Haplotypes A and B (father) and C and D (mother), named like
    /// [`Estimate::haplotypes`]; unknown loci appear as `LOCUS*?`. A and C
    /// are the haplotypes the lowest-numbered typed child received.
    pub haplotypes: [String; 4],
    pub children: Vec<ChildInheritance>,
    /// Loci no inheritance pattern explains (typing error or recombination).
    pub inconsistent_loci: Vec<String>,
    /// Loci where at least one parental allele or its phase is unknown.
    pub ambiguous_loci: Vec<String>,
    pub status: PhaseStatus,
    pub notes: Vec<String>,
    /// Per parent (father, mother), per locus.
    pub(crate) parents: [Vec<ParentLocus>; 2],
}

impl FamilyPhase {
    /// Sample ID to use for each parent: the typed ID, or `family/father (inferred)`.
    pub fn parent_id(&self, p: usize) -> String {
        let (typed, role) = if p == 0 {
            (&self.father, "father")
        } else {
            (&self.mother, "mother")
        };
        typed
            .clone()
            .unwrap_or_else(|| format!("{}/{role} (inferred)", self.family))
    }

    fn failed(fam: &Family, ds: &Dataset, n_loci: usize, reason: &str) -> Self {
        FamilyPhase {
            family: fam.id.clone(),
            father: fam.father.map(|i| ds.samples[i].id.clone()),
            mother: fam.mother.map(|i| ds.samples[i].id.clone()),
            haplotypes: Default::default(),
            children: Vec::new(),
            inconsistent_loci: Vec::new(),
            ambiguous_loci: Vec::new(),
            status: PhaseStatus::Failed(reason.to_string()),
            notes: Vec::new(),
            parents: [
                vec![ParentLocus::Unknown; n_loci],
                vec![ParentLocus::Unknown; n_loci],
            ],
        }
    }
}

/// Can parental alleles x (from the father) and y (from the mother) give child genotype g?
/// An unknown parental allele cannot explain a typed child allele.
fn explains(g: &LocusGeno, x: u32, y: u32) -> bool {
    match (g.a, g.b) {
        (WILD, WILD) => true,
        (a, WILD) => x == a || y == a,
        (a, b) => (x == a && y == b) || (x == b && y == a),
    }
}

/// Ordered (haplotype 1, haplotype 2) allele pairs a parent could have at a locus.
fn parent_pairs(opts: Option<&[(LocusGeno, f64)]>, domain: &[u32]) -> Vec<(u32, u32)> {
    let mut any = domain.to_vec();
    any.push(WILD);
    let typed = opts.filter(|o| o.iter().any(|(g, _)| g.a != WILD));
    let mut pairs = Vec::new();
    match typed {
        None => {
            for &x in &any {
                for &y in &any {
                    pairs.push((x, y));
                }
            }
        }
        Some(opts) => {
            for (g, _) in opts {
                match (g.a, g.b) {
                    (WILD, _) => {}
                    (a, WILD) => {
                        for &d in &any {
                            pairs.push((a, d));
                            pairs.push((d, a));
                        }
                    }
                    (a, b) => {
                        pairs.push((a, b));
                        pairs.push((b, a));
                    }
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

/// Per locus: every (father pair, mother pair) and, per child, which of the
/// four inheritance patterns it allows (bit 2*paternal + maternal).
struct LocusModel {
    quads: Vec<[(u32, u32); 2]>,
    /// Quads grouped by their per-child masks: (masks, quad indices).
    groups: Vec<(Vec<u8>, Vec<usize>)>,
    /// Too many possibilities to search (e.g. an untyped parent and a child
    /// typed at the first field only): the locus fits any inheritance and
    /// is left unphased.
    skipped: bool,
}

/// Most work (parent-pair combinations x child genotype options) spent on one locus.
const MAX_LOCUS_WORK: usize = 50_000_000;

/// How many ordered allele pairs `parent_pairs` would return, without building them.
fn parent_pair_count(opts: Option<&[(LocusGeno, f64)]>, domain: usize) -> usize {
    let any = domain + 1;
    match opts.filter(|o| o.iter().any(|(g, _)| g.a != WILD)) {
        None => any.saturating_mul(any),
        Some(opts) => opts
            .iter()
            .map(|(g, _)| match (g.a, g.b) {
                (WILD, _) => 0,
                (_, WILD) => 2 * any,
                _ => 2,
            })
            .fold(0usize, usize::saturating_add),
    }
}

impl LocusModel {
    fn build(
        father: Option<&[(LocusGeno, f64)]>,
        mother: Option<&[(LocusGeno, f64)]>,
        children: &[&[(LocusGeno, f64)]],
    ) -> Self {
        let mut domain: Vec<u32> = father
            .into_iter()
            .chain(mother)
            .chain(children.iter().copied())
            .flatten()
            .flat_map(|(g, _)| [g.a, g.b])
            .filter(|&c| c != WILD)
            .collect();
        domain.sort_unstable();
        domain.dedup();
        let child_opts: usize = children.iter().map(|c| c.len().max(1)).sum();
        let work = parent_pair_count(father, domain.len())
            .saturating_mul(parent_pair_count(mother, domain.len()))
            .saturating_mul(child_opts);
        if work > MAX_LOCUS_WORK {
            return LocusModel {
                quads: Vec::new(),
                groups: vec![(vec![0b1111; children.len()], Vec::new())],
                skipped: true,
            };
        }
        let (fp, mp) = (parent_pairs(father, &domain), parent_pairs(mother, &domain));
        let mut quads = Vec::with_capacity(fp.len() * mp.len());
        let mut index: FxHashMap<Vec<u8>, usize> = FxHashMap::default();
        let mut groups: Vec<(Vec<u8>, Vec<usize>)> = Vec::new();
        for &f in &fp {
            for &m in &mp {
                let masks: Vec<u8> = children
                    .iter()
                    .map(|opts| {
                        let mut mask = 0u8;
                        for pi in 0..2 {
                            for mi in 0..2 {
                                let x = if pi == 0 { f.0 } else { f.1 };
                                let y = if mi == 0 { m.0 } else { m.1 };
                                if opts.iter().any(|(g, _)| explains(g, x, y)) {
                                    mask |= 1 << (2 * pi + mi);
                                }
                            }
                        }
                        mask
                    })
                    .collect();
                let q = quads.len();
                quads.push([f, m]);
                match index.get(&masks) {
                    Some(&k) => groups[k].1.push(q),
                    None => {
                        index.insert(masks.clone(), groups.len());
                        groups.push((masks, vec![q]));
                    }
                }
            }
        }
        LocusModel {
            quads,
            groups,
            skipped: false,
        }
    }

    /// Does some inheritance pattern fit every child at this locus?
    fn fits(&self) -> bool {
        self.groups
            .iter()
            .any(|(masks, _)| masks.iter().all(|&m| m != 0))
    }

    /// Groups of quads consistent with an inheritance assignment (one 0..4 code per child).
    fn consistent<'a>(&'a self, assign: &'a [u8]) -> impl Iterator<Item = &'a Vec<usize>> + 'a {
        self.groups
            .iter()
            .filter(move |(masks, _)| masks.iter().zip(assign).all(|(&m, &a)| m & (1 << a) != 0))
            .map(|(_, qs)| qs)
    }
}

pub(crate) fn phase_one(
    fam: &Family,
    ds: &Dataset,
    resolved: &[Option<ResolvedSample>],
    codes: &AlleleTable,
) -> FamilyPhase {
    let n_loci = ds.loci.len();
    let typed = |i: Option<usize>| i.and_then(|i| resolved[i].as_ref());
    let parents = [typed(fam.father), typed(fam.mother)];
    let children: Vec<usize> = fam
        .children
        .iter()
        .copied()
        .filter(|&c| resolved[c].is_some())
        .collect();
    if children.is_empty() {
        return FamilyPhase::failed(fam, ds, n_loci, "no typed children");
    }
    if parents.iter().all(Option::is_none) && children.len() < 2 {
        return FamilyPhase::failed(
            fam,
            ds,
            n_loci,
            "no typed parents and fewer than 2 children",
        );
    }
    let mut notes = Vec::new();
    let used = &children[..children.len().min(MAX_PHASING_CHILDREN)];
    if used.len() < children.len() {
        notes.push(format!(
            "phased from the first {} of {} children",
            used.len(),
            children.len()
        ));
    }
    for (p, role) in [(0, "father"), (1, "mother")] {
        if parents[p].is_none() {
            notes.push(format!(
                "{role} not typed: rebuilt from the children where possible"
            ));
        }
    }
    if parents.iter().all(Option::is_none) {
        notes.push("no parents typed: which rebuilt parent is the father is arbitrary".into());
    }

    let models: Vec<LocusModel> = (0..n_loci)
        .map(|l| {
            let kids: Vec<&[(LocusGeno, f64)]> = used
                .iter()
                .map(|&c| resolved[c].as_ref().unwrap().per_locus[l].as_slice())
                .collect();
            LocusModel::build(
                parents[0].map(|s| s.per_locus[l].as_slice()),
                parents[1].map(|s| s.per_locus[l].as_slice()),
                &kids,
            )
        })
        .collect();
    let skipped: Vec<&str> = models
        .iter()
        .zip(&ds.loci)
        .filter(|(m, _)| m.skipped)
        .map(|(_, l)| l.as_str())
        .collect();
    if !skipped.is_empty() {
        notes.push(format!(
            "too many possibilities to phase {} (low-resolution typing with an untyped parent): left unphased",
            skipped.join(", ")
        ));
    }

    // Search inheritance assignments. Children are in label order, so the
    // lowest-numbered typed child defines haplotypes A (father) and C (mother).
    let k = used.len();
    let mut best_score = 0usize;
    let mut best: Vec<Vec<u8>> = Vec::new();
    let mut assign = vec![0u8; k];
    let total = 4usize.pow((k - 1) as u32);
    // With neither parent typed, "father" and "mother" are interchangeable:
    // every pattern has a mirror image with the two swapped. Keep one of each
    // pair, or merging both would make every locus look ambiguous.
    let no_parents = parents.iter().all(Option::is_none);
    let mirror = |a: &[u8]| {
        a.iter()
            .map(|&c| ((c & 1) << 1) | (c >> 1))
            .collect::<Vec<u8>>()
    };
    for n in 0..total {
        let mut rest = n;
        for a in assign.iter_mut().skip(1) {
            *a = (rest % 4) as u8;
            rest /= 4;
        }
        if no_parents && mirror(&assign) < assign {
            continue;
        }
        let score = models
            .iter()
            .filter(|m| m.consistent(&assign).next().is_some())
            .count();
        if score > best_score {
            best_score = score;
            best.clear();
        }
        if score == best_score && best.len() < 256 {
            best.push(assign.clone());
        }
    }

    // Collect what each parent could be at each locus under the best assignments.
    let mut pairs: [Vec<BTreeSet<(u32, u32)>>; 2] =
        [vec![BTreeSet::new(); n_loci], vec![BTreeSet::new(); n_loci]];
    let mut explained = vec![false; n_loci];
    for a in &best {
        for (l, model) in models.iter().enumerate() {
            for qs in model.consistent(a) {
                explained[l] = true;
                for &q in qs {
                    pairs[0][l].insert(model.quads[q][0]);
                    pairs[1][l].insert(model.quads[q][1]);
                }
            }
        }
    }

    // The single value all candidates agree on, or WILD.
    fn unique(mut vals: impl Iterator<Item = u32>) -> u32 {
        let first = vals.next().unwrap_or(WILD);
        if vals.all(|v| v == first) {
            first
        } else {
            WILD
        }
    }
    let mut result: [Vec<ParentLocus>; 2] = [Vec::new(), Vec::new()];
    for p in 0..2 {
        for l in 0..n_loci {
            let typed_here =
                parents[p].is_some_and(|s| s.per_locus[l].iter().any(|(g, _)| g.a != WILD));
            let set = &pairs[p][l];
            let pl = if !explained[l] || set.is_empty() {
                if typed_here {
                    ParentLocus::Unphased
                } else {
                    ParentLocus::Unknown
                }
            } else if set.len() == 1 {
                let &(x, y) = set.iter().next().unwrap();
                if x == WILD && y == WILD {
                    ParentLocus::Unknown
                } else {
                    ParentLocus::Phased(x, y)
                }
            } else if typed_here {
                ParentLocus::Unphased
            } else {
                let x = unique(set.iter().map(|q| q.0));
                let y = unique(set.iter().map(|q| q.1));
                if x == WILD && y == WILD {
                    ParentLocus::Unknown
                } else {
                    ParentLocus::Phased(x, y)
                }
            };
            result[p].push(pl);
        }
    }

    // Name the four haplotypes.
    let allele = |l: usize, c: u32| {
        if c == WILD {
            format!("{}*?", ds.loci[l])
        } else {
            codes.name(AlleleId(c)).to_string()
        }
    };
    let hap = |p: usize, side: usize| -> String {
        (0..n_loci)
            .map(|l| match result[p][l] {
                ParentLocus::Phased(x, y) => allele(l, if side == 0 { x } else { y }),
                _ => allele(l, WILD),
            })
            .collect::<Vec<_>>()
            .join("~")
    };
    let haplotypes = [hap(0, 0), hap(0, 1), hap(1, 0), hap(1, 1)];

    let inconsistent_loci: Vec<String> = (0..n_loci)
        .filter(|&l| !explained[l])
        .map(|l| ds.loci[l].clone())
        .collect();
    let ambiguous_loci: Vec<String> = (0..n_loci)
        .filter(|&l| explained[l])
        .filter(|&l| {
            result.iter().any(|r| match r[l] {
                ParentLocus::Phased(x, y) => x == WILD || y == WILD,
                _ => true,
            })
        })
        .map(|l| ds.loci[l].clone())
        .collect();

    let mut children_out: Vec<ChildInheritance> = children
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            let decided = |shift: u8| -> Option<u8> {
                if i >= k {
                    return None;
                }
                let first = (best[0][i] >> shift) & 1;
                best.iter()
                    .all(|a| (a[i] >> shift) & 1 == first)
                    .then_some(first)
            };
            let [paternal, maternal] = HAPLOTYPE_LABELS;
            ChildInheritance {
                sample_id: ds.samples[c].id.clone(),
                paternal: decided(1).map(|h| paternal[h as usize]),
                maternal: decided(0).map(|h| maternal[h as usize]),
            }
        })
        .collect();

    let status = if inconsistent_loci.is_empty() && ambiguous_loci.is_empty() {
        PhaseStatus::Complete
    } else {
        PhaseStatus::Partial
    };
    if !inconsistent_loci.is_empty() {
        notes.push(format!(
            "no inheritance pattern fits {}: typing error or recombination",
            inconsistent_loci.join(", ")
        ));
    }
    if best.len() > 1 {
        notes.push(format!("{} equally good inheritance patterns", best.len()));
    }
    // A parent who carries the same haplotype twice passes that haplotype to
    // every child, so "which of the two" has no answer; report A (or C).
    for (p, role) in [(0, "father"), (1, "mother")] {
        if !carries_one_haplotype_twice(&result[p]) {
            continue;
        }
        let label = HAPLOTYPE_LABELS[p][0];
        let mut changed = false;
        for child in &mut children_out {
            let side = if p == 0 {
                &mut child.paternal
            } else {
                &mut child.maternal
            };
            if side.is_none() {
                *side = Some(label);
                changed = true;
            }
        }
        if changed {
            let other = HAPLOTYPE_LABELS[p][1];
            notes.push(format!(
                "the {role}'s haplotypes {label} and {other} are identical at every typed locus: \
                 every child is reported as receiving {label}"
            ));
        }
    }
    FamilyPhase {
        family: fam.id.clone(),
        father: fam.father.map(|i| ds.samples[i].id.clone()),
        mother: fam.mother.map(|i| ds.samples[i].id.clone()),
        haplotypes,
        children: children_out,
        inconsistent_loci,
        ambiguous_loci,
        status,
        notes,
        parents: result,
    }
}

/// True when a parent's two haplotypes are known and identical at every locus
/// where anything is known about them (at least one), i.e. the parent carries
/// the same haplotype twice.
fn carries_one_haplotype_twice(loci: &[ParentLocus]) -> bool {
    let mut seen = false;
    for pl in loci {
        match *pl {
            ParentLocus::Phased(x, y) if x != WILD && x == y => seen = true,
            ParentLocus::Unknown => {}
            _ => return false, // different alleles, unknown phase, or one side unknown
        }
    }
    seen
}

/// Put the notes from [`resolve_families`] first on each family's phase.
fn with_fill_notes(mut phases: Vec<FamilyPhase>, notes: Vec<Vec<String>>) -> Vec<FamilyPhase> {
    for (p, mut n) in phases.iter_mut().zip(notes) {
        if !n.is_empty() {
            n.append(&mut p.notes);
            p.notes = n;
        }
    }
    phases
}

/// Phase every family at the target resolution.
pub fn phase_families(
    ds: &Dataset,
    ped: &Pedigree,
    nom: &Nomenclature,
    resolutions: &[Resolution],
) -> Result<Vec<FamilyPhase>, Error> {
    let ((resolver, resolved, _), notes) = resolve_families(ds, ped, nom, resolutions)?;
    Ok(with_fill_notes(
        ped.families
            .iter()
            .map(|f| phase_one(f, ds, &resolved, &resolver.codes))
            .collect(),
        notes,
    ))
}

// ------------------------------------------------------------ estimation

/// Estimate haplotype frequencies from family data, using founders only.
///
/// Each family contributes its two parents (typed or rebuilt from the
/// children); children are left out because their haplotypes are copies of
/// their parents'. With `use_family_phase`, loci the family phased enter the
/// EM with their phase known; otherwise parents enter as unrelated people.
/// A family with no typed parents that cannot be phased contributes its
/// first child. Samples outside any family are used as unrelated people.
///
/// Returns the estimate and the phasing of every family.
pub fn estimate_families(
    ds: &Dataset,
    ped: &Pedigree,
    nom: &Nomenclature,
    resolutions: &[Resolution],
    opts: &EmOptions,
    use_family_phase: bool,
) -> Result<(Estimate, Vec<FamilyPhase>), Error> {
    let clock = std::time::Instant::now();
    let ((resolver, resolved, mut issues), notes) = resolve_families(ds, ped, nom, resolutions)?;
    issues.extend(ped.issues.iter().cloned());
    let phases: Vec<FamilyPhase> = with_fill_notes(
        ped.families
            .iter()
            .map(|f| phase_one(f, ds, &resolved, &resolver.codes))
            .collect(),
        notes,
    );

    let mut samples: Vec<ResolvedSample> = ped
        .unrelated
        .iter()
        .filter_map(|&i| resolved[i].clone())
        .collect();
    for (fam, phase) in ped.families.iter().zip(&phases) {
        let typed = [fam.father, fam.mother].map(|i| i.and_then(|i| resolved[i].as_ref()));
        let before = samples.len();
        for (p, typed_parent) in typed.into_iter().enumerate() {
            if use_family_phase && !matches!(phase.status, PhaseStatus::Failed(_)) {
                let per_locus: Vec<Vec<(LocusGeno, f64)>> = (0..ds.loci.len())
                    .map(|l| match phase.parents[p][l] {
                        ParentLocus::Phased(x, y) => vec![(LocusGeno::phased(x, y), 1.0)],
                        ParentLocus::Unphased => typed_parent.map_or_else(
                            || vec![(LocusGeno::MISSING, 1.0)],
                            |s| s.per_locus[l].clone(),
                        ),
                        ParentLocus::Unknown => vec![(LocusGeno::MISSING, 1.0)],
                    })
                    .collect();
                if per_locus
                    .iter()
                    .any(|o| o.iter().any(|(g, _)| g.a != WILD || g.b != WILD))
                {
                    samples.push(ResolvedSample {
                        id: phase.parent_id(p),
                        per_locus,
                    });
                }
            } else if let Some(s) = typed_parent {
                samples.push(s.clone());
            }
        }
        let first_child = fam.children.iter().find(|&&c| resolved[c].is_some());
        if let (true, Some(&c)) = (samples.len() == before, first_child) {
            samples.push(resolved[c].clone().unwrap());
            issues.push(Issue {
                line: 0,
                sample_id: Some(ds.samples[c].id.clone()),
                message: format!(
                    "family {}: no founders available; this child used instead",
                    fam.id
                ),
            });
        }
    }
    let est = estimate_resolved(
        &ds.loci,
        resolutions,
        &resolver.codes,
        samples,
        issues,
        opts,
        clock,
    )?;
    Ok((est, phases))
}

// ------------------------------------------------------------ comparison

/// One founder: family-phased haplotypes against the EM's most likely pair.
#[derive(Debug, Clone, PartialEq)]
pub struct FounderComparison {
    pub family: String,
    pub sample_id: String,
    pub family_haplotypes: [String; 2],
    pub em_haplotypes: [String; 2],
    pub em_probability: f64,
    /// Both haplotypes identical (in either order).
    pub exact: bool,
    /// Loci where both haplotypes agree, in the better of the two pairings.
    pub loci_agreeing: usize,
    pub n_loci: usize,
}

/// How well an EM estimate matches family phasing.
#[derive(Debug, Clone, PartialEq)]
pub struct FamilyComparison {
    pub founders: Vec<FounderComparison>,
    pub n_exact: usize,
    /// Mean fraction of loci where the EM's best pair agrees with the family.
    pub mean_locus_agreement: f64,
    /// (haplotype, frequency among family-phased founder haplotypes, EM frequency),
    /// most frequent in the families first.
    pub frequencies: Vec<(String, f64, f64)>,
    /// 0.5 * sum |family - EM| over all haplotypes (0 = identical).
    pub total_variation: f64,
    /// Fully phased founder haplotypes counted.
    pub n_family_haplotypes: usize,
}

/// Compare an estimate (usually from [`estimate_families`] with
/// `use_family_phase = false`, i.e. treating founders as unrelated) with the
/// family phasing. Only founders whose two haplotypes are fully known are
/// compared.
pub fn compare_to_families(est: &Estimate, phases: &[FamilyPhase]) -> FamilyComparison {
    let diplotypes: FxHashMap<&str, &crate::estimate::Diplotype> = est
        .diplotypes
        .iter()
        .map(|d| (d.sample_id.as_str(), d))
        .collect();
    let known = |h: &str| !h.is_empty() && !h.contains("*?");
    let mut founders = Vec::new();
    let mut counts: FxHashMap<&str, usize> = FxHashMap::default();
    let mut n_family_haplotypes = 0usize;
    for phase in phases
        .iter()
        .filter(|p| !matches!(p.status, PhaseStatus::Failed(_)))
    {
        for h in &phase.haplotypes {
            if known(h) {
                *counts.entry(h.as_str()).or_default() += 1;
                n_family_haplotypes += 1;
            }
        }
        for p in 0..2 {
            let (h1, h2) = (&phase.haplotypes[2 * p], &phase.haplotypes[2 * p + 1]);
            let id = phase.parent_id(p);
            let (Some(d), true) = (diplotypes.get(id.as_str()), known(h1) && known(h2)) else {
                continue;
            };
            let split = |s: &str| s.split('~').map(str::to_owned).collect::<Vec<_>>();
            let (f1, f2, e1, e2) = (
                split(h1),
                split(h2),
                split(&d.haplotypes[0]),
                split(&d.haplotypes[1]),
            );
            let agree = |a: &[String], b: &[String], c: &[String], e: &[String]| {
                (0..a.len())
                    .filter(|&l| a[l] == c[l] && b[l] == e[l])
                    .count()
            };
            let loci_agreeing = agree(&f1, &f2, &e1, &e2).max(agree(&f1, &f2, &e2, &e1));
            let mut fam_pair = [h1.clone(), h2.clone()];
            fam_pair.sort();
            let mut em_pair = d.haplotypes.clone();
            em_pair.sort();
            founders.push(FounderComparison {
                family: phase.family.clone(),
                sample_id: id,
                exact: fam_pair == em_pair,
                family_haplotypes: [h1.clone(), h2.clone()],
                em_haplotypes: d.haplotypes.clone(),
                em_probability: d.probability,
                loci_agreeing,
                n_loci: f1.len(),
            });
        }
    }
    let n_exact = founders.iter().filter(|f| f.exact).count();
    let mean_locus_agreement = if founders.is_empty() {
        f64::NAN
    } else {
        founders
            .iter()
            .map(|f| f.loci_agreeing as f64 / f.n_loci as f64)
            .sum::<f64>()
            / founders.len() as f64
    };

    let em: FxHashMap<&str, f64> = est
        .haplotypes
        .iter()
        .map(|(h, f)| (h.as_str(), *f))
        .collect();
    let mut keys: Vec<&str> = counts.keys().copied().chain(em.keys().copied()).collect();
    keys.sort_unstable();
    keys.dedup();
    let fam_freq = |h: &str| {
        counts
            .get(h)
            .map_or(0.0, |&c| c as f64 / n_family_haplotypes.max(1) as f64)
    };
    let total_variation = 0.5
        * keys
            .iter()
            .map(|h| (fam_freq(h) - em.get(h).copied().unwrap_or(0.0)).abs())
            .sum::<f64>();
    let mut frequencies: Vec<(String, f64, f64)> = keys
        .iter()
        .map(|h| {
            (
                h.to_string(),
                fam_freq(h),
                em.get(h).copied().unwrap_or(0.0),
            )
        })
        .filter(|t| t.1 > 0.0 || t.2 > 1e-6)
        .collect();
    frequencies.sort_by(|a, b| {
        b.1.total_cmp(&a.1)
            .then(b.2.total_cmp(&a.2))
            .then(a.0.cmp(&b.0))
    });

    FamilyComparison {
        founders,
        n_exact,
        mean_locus_agreement,
        frequencies,
        total_variation,
        n_family_haplotypes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load(table: &str) -> (Dataset, Pedigree) {
        read_family_table_from(table.as_bytes(), &ParseOptions::default()).unwrap()
    }

    fn as_is(ds: &Dataset) -> Vec<Resolution> {
        vec![Resolution::AsItIs; ds.loci.len()]
    }

    const TRIO: &str = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2\n\
        f1,F,01:01,02:01,07:02,08:01\n\
        f1,M,03:01,11:01,35:01,44:02\n\
        f1,C1,01:01,03:01,07:02,35:01\n\
        f1,C2,02:01,11:01,08:01,44:02\n";

    #[test]
    fn roles_and_pedigree() {
        assert_eq!(parse_role("Father"), Some(Role::Father));
        assert_eq!(parse_role("m"), Some(Role::Mother));
        assert_eq!(parse_role("C12"), Some(Role::Child));
        assert_eq!(parse_role("child_2"), Some(Role::Child));
        assert_eq!(parse_role("grandma"), None);
        let (ds, ped) = load(&format!(
            "{TRIO}f1,X,01:01,02:01,07:02,08:01\n,,01:01,01:01,08:01,08:01\n"
        ));
        assert_eq!(ped.families.len(), 1);
        let f = &ped.families[0];
        assert_eq!(
            (f.father, f.mother, f.children.clone()),
            (Some(0), Some(1), vec![2, 3])
        );
        assert_eq!(ped.unrelated, vec![5]);
        assert_eq!(ped.issues.len(), 1, "unknown role X");
        assert_eq!(ds.samples.len(), 6);
    }

    #[test]
    fn trio_is_fully_phased() {
        let (ds, ped) = load(TRIO);
        let phases = phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap();
        let p = &phases[0];
        assert_eq!(p.status, PhaseStatus::Complete);
        assert_eq!(
            p.haplotypes,
            [
                "A*01:01~B*07:02".to_string(),
                "A*02:01~B*08:01".to_string(),
                "A*03:01~B*35:01".to_string(),
                "A*11:01~B*44:02".to_string()
            ]
        );
        assert_eq!(p.children[0].labels(), "AC");
        assert_eq!(p.children[1].paternal, Some('B'));
        assert_eq!(p.children[1].maternal, Some('D'));
    }

    #[test]
    fn lowest_numbered_child_defines_a_and_c_whatever_the_row_order() {
        // Same family as TRIO plus C3, with rows in the order C3, F, C2, M, C1.
        // C1 must get A and C, and the children come back in label order.
        let table = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2\n\
            f1,C3,02:01,03:01,08:01,35:01\n\
            f1,F,01:01,02:01,07:02,08:01\n\
            f1,C2,02:01,11:01,08:01,44:02\n\
            f1,M,03:01,11:01,35:01,44:02\n\
            f1,C1,01:01,03:01,07:02,35:01\n";
        let (ds, ped) = load(table);
        let p = &phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap()[0];
        let order: Vec<&str> = p.children.iter().map(|c| c.sample_id.as_str()).collect();
        assert_eq!(order, ["f1/C1", "f1/C2", "f1/C3"]);
        let labels: Vec<String> = p.children.iter().map(ChildInheritance::labels).collect();
        assert_eq!(labels, ["AC", "BD", "BC"]);
        // A and C are C1's haplotypes.
        assert_eq!(p.haplotypes[0], "A*01:01~B*07:02");
        assert_eq!(p.haplotypes[2], "A*03:01~B*35:01");
        assert_eq!(p.haplotypes[1], "A*02:01~B*08:01");
        assert_eq!(p.haplotypes[3], "A*11:01~B*44:02");
    }

    #[test]
    fn parent_carrying_one_haplotype_twice_gives_a_to_every_child() {
        // The father carries 01:01~07:02 twice; the mother is heterozygous.
        let table = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2\n\
            f1,F,01:01,01:01,07:02,07:02\n\
            f1,M,03:01,11:01,35:01,44:02\n\
            f1,C1,01:01,03:01,07:02,35:01\n\
            f1,C2,01:01,11:01,07:02,44:02\n\
            f1,C3,01:01,03:01,07:02,35:01\n";
        let (ds, ped) = load(table);
        let p = &phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap()[0];
        let labels: Vec<String> = p.children.iter().map(ChildInheritance::labels).collect();
        assert_eq!(labels, ["AC", "AD", "AC"]);
        assert_eq!(p.haplotypes[0], p.haplotypes[1]);
        assert!(
            p.notes
                .iter()
                .any(|n| n.contains("father's haplotypes A and B are identical")),
            "{:?}",
            p.notes
        );
    }

    #[test]
    fn homozygous_copy_that_breaks_mendel_is_undone() {
        // M and C2 each have one A allele typed. Read as homozygous, C2 (02:01
        // twice) cannot come from M (01:01 twice). Undoing either copy fixes it,
        // so both are undone. B shows the copies are not needed elsewhere.
        let table = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2\n\
            f1,F,02:07,02:01,07:02,08:01\n\
            f1,M,01:01,,35:01,44:02\n\
            f1,C1,01:01,02:07,07:02,35:01\n\
            f1,C2,02:01,,08:01,44:02\n\
            f1,C3,01:01,02:01,08:01,35:01\n";
        let (mut ds, ped) = load(table);
        crate::fill_dataset(&mut ds, &crate::FillOptions::default());
        assert_eq!(ds.copied.len(), 2);
        let nom = Nomenclature::default();
        let p = &phase_families(&ds, &ped, &nom, &as_is(&ds)).unwrap()[0];
        assert!(p.inconsistent_loci.is_empty(), "{:?}", p.notes);
        assert!(
            p.notes[0].starts_with("A: blank second allele of f1/M, f1/C2 read as untyped"),
            "{:?}",
            p.notes
        );
        assert!(
            check_mendelian(&ds, &ped, &nom, &as_is(&ds))
                .unwrap()
                .is_empty()
        );

        // When the copy fits the family it stays, with no note.
        let fits = table.replace("f1,C2,02:01,,", "f1,C2,02:01,01:01,");
        let (mut ds, ped) = load(&fits);
        crate::fill_dataset(&mut ds, &crate::FillOptions::default());
        let p = &phase_families(&ds, &ped, &nom, &as_is(&ds)).unwrap()[0];
        assert!(
            !p.notes.iter().any(|n| n.contains("blank second allele")),
            "{:?}",
            p.notes
        );
        assert_eq!(p.haplotypes[2].split('~').next(), Some("A*01:01"));
        assert_eq!(p.haplotypes[3].split('~').next(), Some("A*01:01"));
    }

    #[test]
    fn partly_homozygous_parent_is_resolved_by_other_loci() {
        // The father is homozygous at A (01:01) but not at B: B decides.
        let table = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2\n\
            f1,F,01:01,01:01,07:02,08:01\n\
            f1,M,03:01,11:01,35:01,44:02\n\
            f1,C1,01:01,03:01,07:02,35:01\n\
            f1,C2,01:01,11:01,08:01,44:02\n";
        let (ds, ped) = load(table);
        let p = &phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap()[0];
        let labels: Vec<String> = p.children.iter().map(ChildInheritance::labels).collect();
        assert_eq!(labels, ["AC", "BD"]);
        assert!(!p.notes.iter().any(|n| n.contains("identical")));
    }

    #[test]
    fn child_numbers_sort_numerically() {
        assert_eq!(child_number("C"), 0);
        assert_eq!(child_number("C2"), 2);
        assert_eq!(child_number("C10"), 10);
        assert_eq!(child_number("Child_3"), 3);
        let (_, ped) = load(
            "FAMILY_ID,Family_Member,A_1,A_2\n\
             f1,C10,01:01,01:01\nf1,C2,01:01,01:01\nf1,C,01:01,01:01\nf1,C1,01:01,01:01\n",
        );
        assert_eq!(ped.families[0].children, vec![2, 3, 1, 0]); // C, C1, C2, C10
    }

    #[test]
    fn missing_father_is_rebuilt_from_children() {
        let table = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2\n\
            f1,M,03:01,11:01,35:01,44:02\n\
            f1,C1,01:01,03:01,07:02,35:01\n\
            f1,C2,02:01,11:01,08:01,44:02\n\
            f1,C3,01:01,11:01,07:02,44:02\n";
        let (ds, ped) = load(table);
        let p = &phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap()[0];
        assert_eq!(p.status, PhaseStatus::Complete, "{:?}", p.notes);
        assert_eq!(p.haplotypes[0], "A*01:01~B*07:02");
        assert_eq!(p.haplotypes[1], "A*02:01~B*08:01");
        assert_eq!(p.parent_id(0), "f1/father (inferred)");
    }

    #[test]
    fn untransmitted_haplotype_of_missing_parent_stays_unknown() {
        // One child only: the father's second haplotype is never seen.
        let table = "FAMILY_ID,Family_Member,A_1,A_2\n\
            f1,M,03:01,11:01\n\
            f1,C1,01:01,03:01\n\
            f1,C2,01:01,11:01\n";
        let (ds, ped) = load(table);
        let p = &phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap()[0];
        assert_eq!(p.haplotypes[0], "A*01:01");
        assert_eq!(p.haplotypes[1], "A*?");
        assert_eq!(p.status, PhaseStatus::Partial);
    }

    #[test]
    fn children_only_family_is_phased_up_to_parent_labels() {
        // Father 01:01~07:02 / 02:01~08:01, mother 03:01~35:01 / 11:01~44:02; no parents typed.
        let table = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2\n\
            f1,C1,01:01,03:01,07:02,35:01\n\
            f1,C2,02:01,11:01,08:01,44:02\n\
            f1,C3,01:01,11:01,07:02,44:02\n";
        let (ds, ped) = load(table);
        let p = &phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap()[0];
        assert_eq!(p.status, PhaseStatus::Complete, "{:?}", p.notes);
        let mut parents = [
            [p.haplotypes[0].clone(), p.haplotypes[1].clone()],
            [p.haplotypes[2].clone(), p.haplotypes[3].clone()],
        ];
        parents.iter_mut().for_each(|h| h.sort());
        parents.sort();
        assert_eq!(
            parents,
            [
                ["A*01:01~B*07:02".to_string(), "A*02:01~B*08:01".to_string()],
                ["A*03:01~B*35:01".to_string(), "A*11:01~B*44:02".to_string()]
            ]
        );
    }

    #[test]
    fn mendel_errors_are_reported() {
        let table = format!("{TRIO}f1,C3,24:02,03:01,07:02,35:01\n");
        let (ds, ped) = load(&table);
        let issues = check_mendelian(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap();
        assert_eq!(issues.len(), 1);
        assert_eq!(issues[0].sample_id.as_deref(), Some("f1/C3"));
        assert_eq!(issues[0].locus.as_deref(), Some("A"));
        // The phasing still works and flags nothing for the other loci.
        let p = &phase_families(&ds, &ped, &Nomenclature::default(), &as_is(&ds)).unwrap()[0];
        assert_eq!(p.inconsistent_loci, vec!["A".to_string()]);
    }

    #[test]
    fn founders_only_estimate_with_phase_counts_haplotypes() {
        // Two trio families: founders carry 8 haplotypes, all phased.
        let table = format!(
            "{TRIO}f2,F,01:01,02:01,07:02,08:01\nf2,M,01:01,24:02,07:02,15:01\n\
             f2,C1,01:01,01:01,07:02,07:02\nf2,C2,02:01,24:02,08:01,15:01\n"
        );
        let (ds, ped) = load(&table);
        let opts = EmOptions {
            epsilon: 1e-9,
            ..Default::default()
        };
        let (est, phases) = estimate_families(
            &ds,
            &ped,
            &Nomenclature::default(),
            &as_is(&ds),
            &opts,
            true,
        )
        .unwrap();
        assert!(phases.iter().all(|p| p.status == PhaseStatus::Complete));
        assert_eq!(est.n_samples, 4, "founders only");
        let f = |h: &str| {
            est.haplotypes
                .iter()
                .find(|x| x.0 == h)
                .map_or(0.0, |x| x.1)
        };
        assert!((f("A*01:01~B*07:02") - 3.0 / 8.0).abs() < 1e-9);
        assert!((f("A*02:01~B*08:01") - 2.0 / 8.0).abs() < 1e-9);
        assert!((f("A*24:02~B*15:01") - 1.0 / 8.0).abs() < 1e-9);

        // Unphased founders give a comparable estimate; compare it to the families.
        let (unphased, _) = estimate_families(
            &ds,
            &ped,
            &Nomenclature::default(),
            &as_is(&ds),
            &opts,
            false,
        )
        .unwrap();
        let cmp = compare_to_families(&unphased, &phases);
        assert_eq!(cmp.founders.len(), 4);
        assert_eq!(cmp.n_family_haplotypes, 8);
        assert!(cmp.mean_locus_agreement > 0.0);
        let phased_cmp = compare_to_families(&est, &phases);
        assert_eq!(
            phased_cmp.n_exact, 4,
            "phased input reproduces the families"
        );
        assert!(phased_cmp.total_variation < 1e-9);
    }
}
