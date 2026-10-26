//! Pedigree checks: find members whose typing does not fit their stated
//! relationships.
//!
//! [`audit_families`] runs five checks on every family, at the target
//! resolution, and returns one [`PedigreeFlag`] per problem:
//!
//! 1. **Parent exclusion.** A child that shares no allele with a parent at two
//!    or more loci is not that parent's child (non-paternity, sample swap, wrong
//!    label). A single such locus is reported separately: it is more often a
//!    typing error, a null allele or a mutation.
//! 2. **Trio consistency.** Each parent fits on its own, but no combination of
//!    one allele from each explains the child.
//! 3. **Duplicates.** A child identical to a parent, or two parents identical to
//!    each other, at almost every locus: usually a duplicated row. Siblings are
//!    not compared (HLA-identical siblings are expected). Identical people in
//!    different families are also reported.
//! 4. **Siblings.** With a parent untyped, the children together may carry more
//!    alleles than the parents could pass on. The child whose removal fixes
//!    this is reported as not a full sibling.
//! 5. **Haplotype fit.** Each child is compared with the parental haplotypes
//!    phased from the *rest* of the family. One switch between a parent's two
//!    haplotypes along chromosome 6 is reported as a probable recombination
//!    (with the interval); more switches, or alleles a rebuilt parent cannot
//!    have, point to typing errors or a wrong relationship.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::Error;
use crate::dataset::Dataset;
use crate::estimate::{LocusGeno, ResolvedSample, WILD, resolve_dataset};
use crate::family::{
    Family, ParentLocus, Pedigree, PhaseStatus, child_possible, options_text, parent_has, phase_one,
};
use crate::intern::{AlleleId, AlleleTable};
use crate::nomenclature::Nomenclature;
use crate::resolve::Resolution;

/// Loci in chromosome 6 order (telomere to centromere), for recombination checks.
pub const GENOMIC_ORDER: &[&str] = &[
    "HFE", "F", "V", "G", "H", "K", "A", "J", "L", "E", "C", "B", "MICA", "MICB", "DRA", "DRB9",
    "DRB5", "DRB4", "DRB3", "DRB1", "DQA1", "DQB1", "DQA2", "DQB2", "DOB", "TAP2", "TAP1", "DMB",
    "DMA", "DOA", "DPA1", "DPB1",
];

/// Kind of pedigree problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum FlagKind {
    /// Family structure: unknown role, second parent, unreadable typing.
    Structure,
    /// Child shares no allele with a parent at several loci.
    ParentExcluded,
    /// One locus does not fit (typing error, null allele, mutation).
    SingleLocusMismatch,
    /// Each parent fits alone, but not both together.
    TrioInconsistent,
    /// Near-identical to a parent or spouse (or identical person in another family).
    PossibleDuplicate,
    /// Child carries alleles the other children's parents cannot also have passed on.
    NotFullSibling,
    /// Child's inheritance from one parent switches once along chromosome 6.
    Recombination,
    /// Child does not fit the family's haplotypes at several loci.
    InheritanceMismatch,
}

impl FlagKind {
    pub fn as_str(self) -> &'static str {
        match self {
            FlagKind::Structure => "structure",
            FlagKind::ParentExcluded => "parent_excluded",
            FlagKind::SingleLocusMismatch => "single_locus_mismatch",
            FlagKind::TrioInconsistent => "trio_inconsistent",
            FlagKind::PossibleDuplicate => "possible_duplicate",
            FlagKind::NotFullSibling => "not_full_sibling",
            FlagKind::Recombination => "recombination",
            FlagKind::InheritanceMismatch => "inheritance_mismatch",
        }
    }
}

/// One problem found in a family.
#[derive(Debug, Clone, PartialEq)]
pub struct PedigreeFlag {
    pub family: String,
    /// The member the problem is about (usually the child).
    pub sample_id: Option<String>,
    /// The other member involved (e.g. the excluded parent), if any.
    pub relative: Option<String>,
    pub kind: FlagKind,
    pub loci: Vec<String>,
    pub message: String,
}

/// Thresholds for [`audit_families`].
#[derive(Debug, Clone)]
pub struct AuditOptions {
    /// Loci without a shared allele needed to call a parent excluded (default 2).
    pub min_exclusions: usize,
    /// Fraction of identical loci at which two members look duplicated (default 0.9).
    pub duplicate_fraction: f64,
    /// Fewest loci compared before calling a duplicate (default 5).
    pub min_loci_compared: usize,
    /// Also report people in different families with identical typing at every
    /// locus (default off: unrelated people with common haplotypes often match).
    pub cross_family_duplicates: bool,
}

impl Default for AuditOptions {
    fn default() -> Self {
        AuditOptions {
            min_exclusions: 2,
            duplicate_fraction: 0.9,
            min_loci_compared: 5,
            cross_family_duplicates: false,
        }
    }
}

type Opts<'a> = &'a [(LocusGeno, f64)];

/// Both alleles known in every possible genotype.
fn fully_typed(o: Opts) -> bool {
    !o.is_empty() && o.iter().all(|(g, _)| g.a != WILD && g.b != WILD)
}

/// The child shares no allele with the parent (both fully typed).
fn excluded(child: Opts, parent: Opts) -> bool {
    fully_typed(child)
        && fully_typed(parent)
        && !child.iter().any(|(c, _)| {
            parent
                .iter()
                .any(|(p, _)| parent_has(p, c.a) || parent_has(p, c.b))
        })
}

fn trio_possible(child: Opts, father: Opts, mother: Opts) -> bool {
    child.iter().any(|(c, _)| {
        father
            .iter()
            .any(|(f, _)| mother.iter().any(|(m, _)| child_possible(c, f, m)))
    })
}

fn same_genotype(x: Opts, y: Opts) -> bool {
    x.iter().any(|(a, _)| y.iter().any(|(b, _)| a == b))
}

/// Could these children all be full siblings at one locus? Each child passes one
/// allele from each parent; a typed parent limits what it can pass, an untyped
/// parent can pass at most two different alleles.
fn siblings_possible(kids: &[Opts], father: Option<Opts>, mother: Option<Opts>) -> bool {
    let allowed = |p: Option<Opts>| -> Option<FxHashSet<u32>> {
        let p = p.filter(|o| o.iter().any(|(g, _)| g.a != WILD))?;
        if p.iter().any(|(g, _)| g.b == WILD) {
            return None; // half-typed parent: anything could fit
        }
        Some(p.iter().flat_map(|(g, _)| [g.a, g.b]).collect())
    };
    let (fa, ma) = (allowed(father), allowed(mother));
    let kids: Vec<Opts> = kids
        .iter()
        .copied()
        .filter(|k| k.iter().any(|(g, _)| g.a != WILD))
        .collect();

    fn search(
        kids: &[Opts],
        i: usize,
        fs: &mut Vec<u32>,
        ms: &mut Vec<u32>,
        fa: &Option<FxHashSet<u32>>,
        ma: &Option<FxHashSet<u32>>,
    ) -> bool {
        let Some(opts) = kids.get(i) else { return true };
        let fits = |x: u32, set: &Vec<u32>, allowed: &Option<FxHashSet<u32>>| {
            x == WILD
                || match allowed {
                    Some(a) => a.contains(&x),
                    None => set.contains(&x) || set.len() < 2,
                }
        };
        for (g, _) in opts.iter() {
            for (x, y) in [(g.a, g.b), (g.b, g.a)] {
                if fits(x, fs, fa) && fits(y, ms, ma) {
                    let (fnew, mnew) =
                        (x != WILD && !fs.contains(&x), y != WILD && !ms.contains(&y));
                    if fnew {
                        fs.push(x);
                    }
                    if mnew {
                        ms.push(y);
                    }
                    let ok = search(kids, i + 1, fs, ms, fa, ma);
                    if fnew {
                        fs.pop();
                    }
                    if mnew {
                        ms.pop();
                    }
                    if ok {
                        return true;
                    }
                }
            }
        }
        false
    }
    search(&kids, 0, &mut Vec::new(), &mut Vec::new(), &fa, &ma)
}

/// Loci in chromosome order (unknown loci last, in input order).
fn genomic_order(loci: &[String]) -> Vec<usize> {
    let pos = |l: &str| {
        GENOMIC_ORDER
            .iter()
            .position(|g| *g == l)
            .unwrap_or(usize::MAX)
    };
    let mut order: Vec<usize> = (0..loci.len()).collect();
    order.sort_by_key(|&i| (pos(&loci[i]), i));
    order
}

/// Check every family for relationships the typing does not support.
pub fn audit_families(
    ds: &Dataset,
    ped: &Pedigree,
    nom: &Nomenclature,
    resolutions: &[Resolution],
    opts: &AuditOptions,
) -> Result<Vec<PedigreeFlag>, Error> {
    let (resolver, resolved, issues) = resolve_dataset(ds, nom, resolutions)?;
    let codes = &resolver.codes;
    let mut family_of: FxHashMap<usize, &str> = FxHashMap::default();
    for f in &ped.families {
        for i in f.father.iter().chain(&f.mother).chain(&f.children) {
            family_of.insert(*i, &f.id);
        }
    }
    let mut flags = Vec::new();
    for issue in ped.issues.iter() {
        flags.push(PedigreeFlag {
            family: String::new(),
            sample_id: issue.sample_id.clone(),
            relative: None,
            kind: FlagKind::Structure,
            loci: Vec::new(),
            message: issue.message.clone(),
        });
    }
    for issue in &issues {
        let fam = issue
            .sample_id
            .as_ref()
            .and_then(|s| ds.samples.iter().position(|x| &x.id == s))
            .and_then(|i| family_of.get(&i).copied())
            .unwrap_or_default();
        flags.push(PedigreeFlag {
            family: fam.to_string(),
            sample_id: issue.sample_id.clone(),
            relative: None,
            kind: FlagKind::Structure,
            loci: Vec::new(),
            message: format!("typing not usable: {}", issue.message),
        });
    }
    let order = genomic_order(&ds.loci);
    for fam in &ped.families {
        audit_family(fam, ds, &resolved, codes, opts, &order, &mut flags);
    }

    // Exact duplicates across families: the same typing at every locus.
    let mut seen: FxHashMap<Vec<LocusGeno>, usize> = FxHashMap::default();
    let candidates = if opts.cross_family_duplicates {
        resolved.len()
    } else {
        0
    };
    for (i, s) in resolved.iter().enumerate().take(candidates) {
        let Some(s) = s else { continue };
        let typed = s.per_locus.iter().filter(|o| fully_typed(o)).count();
        if typed < opts.min_loci_compared || s.per_locus.iter().any(|o| o.len() != 1) {
            continue;
        }
        let key: Vec<LocusGeno> = s.per_locus.iter().map(|o| o[0].0).collect();
        match seen.get(&key) {
            Some(&j) if family_of.get(&j) != family_of.get(&i) => flags.push(PedigreeFlag {
                family: family_of.get(&i).copied().unwrap_or_default().to_string(),
                sample_id: Some(ds.samples[i].id.clone()),
                relative: Some(ds.samples[j].id.clone()),
                kind: FlagKind::PossibleDuplicate,
                loci: Vec::new(),
                message: format!(
                    "identical at all {typed} typed loci to {} in another family: the same person, a duplicated \
                     sample, or (with common haplotypes) two unrelated people who match by chance",
                    ds.samples[j].id
                ),
            }),
            Some(_) => {}
            None => {
                seen.insert(key, i);
            }
        }
    }
    Ok(flags)
}

#[allow(clippy::too_many_arguments)]
fn audit_family(
    fam: &Family,
    ds: &Dataset,
    resolved: &[Option<ResolvedSample>],
    codes: &AlleleTable,
    opts: &AuditOptions,
    order: &[usize],
    flags: &mut Vec<PedigreeFlag>,
) {
    let n_loci = ds.loci.len();
    let id = |i: usize| ds.samples[i].id.clone();
    let names = |ls: &[usize]| ls.iter().map(|&l| ds.loci[l].clone()).collect::<Vec<_>>();
    let mut flag = |sample: usize,
                    relative: Option<usize>,
                    kind: FlagKind,
                    loci: &[usize],
                    message: String| {
        flags.push(PedigreeFlag {
            family: fam.id.clone(),
            sample_id: Some(id(sample)),
            relative: relative.map(id),
            kind,
            loci: names(loci),
            message,
        })
    };
    let typed = |i: Option<usize>| i.and_then(|i| resolved[i].as_ref().map(|s| (i, s)));
    let parents = [typed(fam.father), typed(fam.mother)];
    let role = ["father", "mother"];
    let kids: Vec<(usize, &ResolvedSample)> = fam
        .children
        .iter()
        .filter_map(|&c| resolved[c].as_ref().map(|s| (c, s)))
        .collect();
    // Children already explained by an exclusion, per parent.
    let mut excluded_from: FxHashSet<(usize, usize)> = FxHashSet::default();
    // Flags not about one member; added after the per-member ones.
    let mut family_level: Vec<PedigreeFlag> = Vec::new();
    // Children named together in an ambiguous sibling flag (not flagged one by one).
    let mut ambiguous_kids: FxHashSet<usize> = FxHashSet::default();
    // Children whose row looks copied from a parent.
    let mut duplicate_rows: Vec<usize> = Vec::new();
    let mut flagged_loci: FxHashMap<usize, FxHashSet<usize>> = FxHashMap::default();

    // 1-2. Each child against each parent, then against both together.
    for &(c, cs) in &kids {
        for p in 0..2 {
            let Some((pi, ps)) = parents[p] else { continue };
            let loci: Vec<usize> = (0..n_loci)
                .filter(|&l| excluded(&cs.per_locus[l], &ps.per_locus[l]))
                .collect();
            let compared = (0..n_loci)
                .filter(|&l| fully_typed(&cs.per_locus[l]) && fully_typed(&ps.per_locus[l]))
                .count();
            let detail = || {
                loci.iter()
                    .map(|&l| {
                        format!(
                            "{} child {} / {} {}",
                            ds.loci[l],
                            options_text(&cs.per_locus[l], codes),
                            role[p],
                            options_text(&ps.per_locus[l], codes)
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("; ")
            };
            if loci.len() >= opts.min_exclusions {
                excluded_from.insert((c, p));
                flag(
                    c,
                    Some(pi),
                    FlagKind::ParentExcluded,
                    &loci,
                    format!(
                        "shares no allele with the {} at {} of {compared} loci compared: not the biological {}, or a sample swap ({})",
                        role[p],
                        loci.len(),
                        role[p],
                        detail()
                    ),
                );
            } else if loci.len() == 1 {
                flag(
                    c,
                    Some(pi),
                    FlagKind::SingleLocusMismatch,
                    &loci,
                    format!(
                        "shares no allele with the {} at one locus: typing error, null allele or mutation ({})",
                        role[p],
                        detail()
                    ),
                );
            }
            flagged_loci.entry(c).or_default().extend(loci);
        }
        if let (Some((_, f)), Some((_, m))) = (parents[0], parents[1]) {
            let done = flagged_loci.entry(c).or_default();
            let loci: Vec<usize> = (0..n_loci)
                .filter(|l| !done.contains(l))
                .filter(|&l| !trio_possible(&cs.per_locus[l], &f.per_locus[l], &m.per_locus[l]))
                .collect();
            if !loci.is_empty() {
                done.extend(loci.iter().copied());
                let why = if loci.len() >= opts.min_exclusions {
                    "wrong parents or a sample swap"
                } else {
                    "typing error, null allele or mutation"
                };
                let detail: Vec<String> = loci
                    .iter()
                    .map(|&l| {
                        format!(
                            "{} child {} / father {} / mother {}",
                            ds.loci[l],
                            options_text(&cs.per_locus[l], codes),
                            options_text(&f.per_locus[l], codes),
                            options_text(&m.per_locus[l], codes)
                        )
                    })
                    .collect();
                flag(
                    c,
                    None,
                    FlagKind::TrioInconsistent,
                    &loci,
                    format!(
                        "each parent fits alone but not both together at {} loci: {why} ({})",
                        loci.len(),
                        detail.join("; ")
                    ),
                );
            }
        }
    }

    // 3. Near-identical parent-child or parent-parent pairs.
    let mut members: Vec<(usize, &ResolvedSample, bool)> = parents
        .iter()
        .flatten()
        .map(|&(i, s)| (i, s, true))
        .collect();
    members.extend(kids.iter().map(|&(i, s)| (i, s, false)));
    for a in 0..members.len() {
        for b in a + 1..members.len() {
            let ((i, si, parent_i), (j, sj, parent_j)) = (members[a], members[b]);
            if !parent_i && !parent_j {
                continue; // identical siblings are expected
            }
            let compared: Vec<usize> = (0..n_loci)
                .filter(|&l| fully_typed(&si.per_locus[l]) && fully_typed(&sj.per_locus[l]))
                .collect();
            let same = compared
                .iter()
                .filter(|&&l| same_genotype(&si.per_locus[l], &sj.per_locus[l]))
                .count();
            if compared.len() >= opts.min_loci_compared
                && same as f64 / compared.len() as f64 >= opts.duplicate_fraction
            {
                // A child can match a parent naturally (the other parent passed on a
                // haplotype identical to this parent's second one). A copied row
                // almost never fits the other parent too, so with both parents
                // typed only a match that fails the trio check is reported.
                let both_parents_typed = parents.iter().all(Option::is_some);
                let other_parent_misfit = flagged_loci.get(&j).filter(|s| !s.is_empty());
                let why = match (parent_i && parent_j, both_parents_typed, other_parent_misfit) {
                    (true, _, _) => "the two parents look like the same person: probably a duplicated row".to_string(),
                    (false, true, None) => continue, // natural parent-child match
                    (false, true, Some(bad)) => {
                        let mut bad: Vec<usize> = bad.iter().copied().collect();
                        bad.sort_unstable();
                        format!(
                            "and the child does not fit both parents together at {}: probably the parent's row copied into the child's",
                            names(&bad).join(", ")
                        )
                    }
                    (false, false, _) => "probably a duplicated row; the other parent is untyped, so a natural match \
                                          (both parents sharing a haplotype) cannot be ruled out"
                        .to_string(),
                };
                let differ: Vec<usize> = compared
                    .iter()
                    .copied()
                    .filter(|&l| !same_genotype(&si.per_locus[l], &sj.per_locus[l]))
                    .collect();
                if !parent_j {
                    duplicate_rows.push(j);
                }
                flag(
                    j,
                    Some(i),
                    FlagKind::PossibleDuplicate,
                    &differ,
                    format!(
                        "identical to {} at {same} of {} loci{}, {why}",
                        id(i),
                        compared.len(),
                        if differ.is_empty() {
                            String::new()
                        } else {
                            format!(" (differs only at {})", names(&differ).join(", "))
                        }
                    ),
                );
            }
        }
    }

    // 4. Siblings: with a parent untyped, can the children share two parents?
    if kids.len() >= 2 && parents.iter().any(Option::is_none) {
        let mut bad: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
        let mut unattributed: Vec<usize> = Vec::new();
        let mut candidates: Vec<usize> = Vec::new();
        for l in 0..n_loci {
            let po = |p: usize| parents[p].map(|(_, s)| s.per_locus[l].as_slice());
            let all: Vec<Opts> = kids
                .iter()
                .map(|(_, s)| s.per_locus[l].as_slice())
                .collect();
            if siblings_possible(&all, po(0), po(1)) {
                continue;
            }
            let culprits: Vec<usize> = (0..kids.len())
                .filter(|&k| {
                    let rest: Vec<Opts> = all
                        .iter()
                        .enumerate()
                        .filter(|&(i, _)| i != k)
                        .map(|(_, o)| *o)
                        .collect();
                    siblings_possible(&rest, po(0), po(1))
                })
                .collect();
            match culprits.as_slice() {
                [k] => bad.entry(kids[*k].0).or_default().push(l),
                many => {
                    unattributed.push(l);
                    candidates.extend(many.iter().map(|&k| kids[k].0));
                }
            }
        }
        let mut sorted: Vec<_> = bad.into_iter().collect();
        sorted.sort();
        for (c, loci) in sorted {
            if excluded_from.iter().any(|&(x, _)| x == c) {
                continue; // already explained
            }
            let kind = if loci.len() >= opts.min_exclusions {
                FlagKind::NotFullSibling
            } else {
                FlagKind::SingleLocusMismatch
            };
            flag(
                c,
                None,
                kind,
                &loci,
                format!(
                    "carries alleles at {} that the other children's parents could not also have passed on: {}",
                    names(&loci).join(", "),
                    if kind == FlagKind::NotFullSibling {
                        "half-sibling, unrelated child or sample swap"
                    } else {
                        "typing error, or not a full sibling"
                    }
                ),
            );
        }
        candidates.sort_unstable();
        candidates.dedup();
        ambiguous_kids.extend(candidates.iter().copied());
        if !unattributed.is_empty() {
            let who = if candidates.is_empty() {
                "and no single child explains it".to_string()
            } else {
                format!(
                    "and removing any one of {} would fix it: one of them is probably a half-sibling, \
                     but without the untyped parent the data cannot say which",
                    candidates
                        .iter()
                        .map(|&c| id(c))
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            };
            family_level.push(PedigreeFlag {
                family: fam.id.clone(),
                sample_id: None,
                relative: None,
                kind: FlagKind::NotFullSibling,
                loci: names(&unattributed),
                message: format!(
                    "the children cannot all be full siblings at {}, {who}",
                    names(&unattributed).join(", ")
                ),
            });
        }
    }

    // 5. Each child against haplotypes phased from the rest of the family.
    // Children already shown not to belong (excluded parent, copied row) are
    // left out of the reference phasing, so they cannot distort it.
    let unreliable: FxHashSet<usize> = excluded_from
        .iter()
        .map(|&(c, _)| c)
        .chain(duplicate_rows.iter().copied())
        .collect();
    for &(c, cs) in &kids {
        let others: Vec<usize> = fam
            .children
            .iter()
            .copied()
            .filter(|&x| x != c && !unreliable.contains(&x))
            .collect();
        let n_others = others.iter().filter(|&&x| resolved[x].is_some()).count();
        let n_parents = parents.iter().flatten().count();
        let informative = match n_parents {
            2 => n_others >= 1,
            1 => n_others >= 2,
            _ => n_others >= 3,
        };
        if !informative {
            continue;
        }
        let rest = Family {
            id: fam.id.clone(),
            father: fam.father,
            mother: fam.mother,
            children: others,
        };
        let phase = phase_one(&rest, ds, resolved, codes);
        if matches!(phase.status, PhaseStatus::Failed(_)) {
            continue;
        }
        for p in 0..2 {
            if excluded_from.contains(&(c, p)) {
                continue;
            }
            let rebuilt = parents[p].is_none();
            let mut seq: Vec<(usize, u8)> = Vec::new();
            let mut nonfit: Vec<usize> = Vec::new();
            for &l in order {
                let ParentLocus::Phased(x, y) = phase.parents[p][l] else {
                    continue;
                };
                if x == WILD || y == WILD || x == y {
                    continue;
                }
                let [(g, _)] = cs.per_locus[l].as_slice() else {
                    continue;
                };
                if g.a == WILD {
                    continue;
                }
                let has = |v: u32| g.a == v || g.b == v;
                match (has(x), has(y)) {
                    (true, false) => seq.push((l, 0)),
                    (false, true) => seq.push((l, 1)),
                    (false, false) if rebuilt && g.b != WILD => nonfit.push(l),
                    _ => {}
                }
            }
            let switches: Vec<usize> = (1..seq.len())
                .filter(|&i| seq[i].1 != seq[i - 1].1)
                .collect();
            let label = |s: &[(usize, u8)]| {
                s.iter()
                    .map(|&(l, _)| ds.loci[l].as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            match switches.as_slice() {
                [] => {}
                [k] => {
                    let (before, after) = seq.split_at(*k);
                    let short = before.len() == 1 || after.len() == 1;
                    flag(
                        c,
                        parents[p].map(|(i, _)| i),
                        FlagKind::Recombination,
                        &[before.last().unwrap().0, after[0].0],
                        format!(
                            "inherits the {}'s haplotype {} at {} and haplotype {} at {}: crossover between {} and {}{}",
                            role[p],
                            before[0].1 + 1,
                            label(before),
                            after[0].1 + 1,
                            label(after),
                            ds.loci[before.last().unwrap().0],
                            ds.loci[after[0].0],
                            if short {
                                " (or a typing error at the single locus on one side)"
                            } else {
                                ""
                            }
                        ),
                    );
                }
                many => flag(
                    c,
                    parents[p].map(|(i, _)| i),
                    FlagKind::InheritanceMismatch,
                    &many.iter().map(|&i| seq[i].0).collect::<Vec<_>>(),
                    format!(
                        "inheritance from the {} switches {} times along chromosome 6 ({}): typing errors, or not this {}'s child",
                        role[p],
                        many.len(),
                        label(&seq),
                        role[p]
                    ),
                ),
            }
            let fresh: Vec<usize> = if ambiguous_kids.contains(&c) {
                Vec::new() // already reported in the family-level sibling flag
            } else {
                nonfit
                    .into_iter()
                    .filter(|l| !flagged_loci.get(&c).is_some_and(|s| s.contains(l)))
                    .collect()
            };
            if !fresh.is_empty() {
                let kind = if fresh.len() >= opts.min_exclusions {
                    FlagKind::InheritanceMismatch
                } else {
                    FlagKind::SingleLocusMismatch
                };
                let allele_pair = |l: usize| match phase.parents[p][l] {
                    ParentLocus::Phased(x, y) => {
                        format!("{}/{}", codes.name(AlleleId(x)), codes.name(AlleleId(y)))
                    }
                    _ => "?".into(),
                };
                flag(
                    c,
                    None,
                    kind,
                    &fresh,
                    format!(
                        "carries neither allele the untyped {} must have (from the other children) at {}: {}",
                        role[p],
                        fresh
                            .iter()
                            .map(|&l| format!("{} {}", ds.loci[l], allele_pair(l)))
                            .collect::<Vec<_>>()
                            .join(", "),
                        if kind == FlagKind::InheritanceMismatch {
                            "half-sibling, unrelated child or sample swap"
                        } else {
                            "typing error, or not a full sibling"
                        }
                    ),
                );
            }
        }
    }
    flags.extend(family_level);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::family::read_family_table_from;
    use crate::parse::ParseOptions;

    fn audit(table: &str) -> Vec<PedigreeFlag> {
        let (ds, ped) = read_family_table_from(table.as_bytes(), &ParseOptions::default()).unwrap();
        let res = vec![Resolution::AsItIs; ds.loci.len()];
        // The test tables have 3 loci; allow duplicate calls on that few.
        let opts = AuditOptions {
            min_loci_compared: 3,
            cross_family_duplicates: true,
            ..AuditOptions::default()
        };
        audit_families(&ds, &ped, &Nomenclature::default(), &res, &opts).unwrap()
    }

    fn kinds(flags: &[PedigreeFlag], sample: &str) -> Vec<FlagKind> {
        flags
            .iter()
            .filter(|f| f.sample_id.as_deref() == Some(sample))
            .map(|f| f.kind)
            .collect()
    }

    // Father 01~07~03 / 02~08~04, mother 03~35~11 / 11~44~13 (A~B~DRB1).
    const HEAD: &str = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2,DRB1_1,DRB1_2\n\
        f1,F,01:01,02:01,07:02,08:01,03:01,04:01\n\
        f1,M,03:01,11:01,35:01,44:02,11:01,13:01\n\
        f1,C1,01:01,03:01,07:02,35:01,03:01,11:01\n\
        f1,C2,02:01,11:01,08:01,44:02,04:01,13:01\n";

    #[test]
    fn clean_family_has_no_flags() {
        assert!(audit(HEAD).is_empty());
    }

    #[test]
    fn non_paternity_is_flagged() {
        // C3 has the mother's haplotype but another father.
        let flags = audit(&format!(
            "{HEAD}f1,C3,24:02,03:01,15:01,35:01,15:01,11:01\n"
        ));
        assert_eq!(kinds(&flags, "f1/C3"), [FlagKind::ParentExcluded]);
        assert_eq!(flags[0].relative.as_deref(), Some("f1/F"));
        assert_eq!(flags[0].loci, ["A", "B", "DRB1"]);
    }

    #[test]
    fn single_locus_mismatch_is_not_called_exclusion() {
        let flags = audit(&format!(
            "{HEAD}f1,C3,24:02,03:01,07:02,35:01,03:01,11:01\n"
        ));
        assert_eq!(kinds(&flags, "f1/C3"), [FlagKind::SingleLocusMismatch]);
    }

    #[test]
    fn child_duplicating_a_parent_is_flagged() {
        let flags = audit(&format!(
            "{HEAD}f1,C3,03:01,11:01,35:01,44:02,11:01,13:01\n"
        ));
        let k = kinds(&flags, "f1/C3");
        assert!(k.contains(&FlagKind::PossibleDuplicate), "{flags:?}");
    }

    #[test]
    fn recombinant_child_is_flagged_with_the_interval() {
        // C3: from the father 01:01~07:02 (hap 1) at A, B but 04:01 (hap 2) at DRB1.
        let flags = audit(&format!(
            "{HEAD}f1,C3,01:01,03:01,07:02,35:01,04:01,11:01\n"
        ));
        let rec: Vec<_> = flags
            .iter()
            .filter(|f| f.kind == FlagKind::Recombination)
            .collect();
        assert_eq!(rec.len(), 1, "{flags:?}");
        assert_eq!(rec[0].sample_id.as_deref(), Some("f1/C3"));
        assert_eq!(rec[0].loci, ["B", "DRB1"]);
        assert!(!flags.iter().any(|f| f.kind == FlagKind::ParentExcluded));
    }

    const NO_FATHER: &str = "FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2,DRB1_1,DRB1_2\n\
        f1,M,03:01,11:01,35:01,44:02,11:01,13:01\n\
        f1,C1,01:01,03:01,07:02,35:01,03:01,11:01\n\
        f1,C2,02:01,11:01,08:01,44:02,04:01,13:01\n\
        f1,C4,01:01,11:01,07:02,44:02,03:01,13:01\n";
    const HALF_SIB: &str = "f1,C3,24:02,03:01,15:01,35:01,15:01,11:01\n";

    #[test]
    fn half_sibling_found_without_the_father() {
        // Father untyped; both his haplotypes are seen twice (C5 repeats C2's),
        // so only C3, with a third paternal haplotype, can be the half-sibling.
        let table = format!("{NO_FATHER}f1,C5,02:01,03:01,08:01,35:01,04:01,11:01\n{HALF_SIB}");
        let flags = audit(&table);
        let k = kinds(&flags, "f1/C3");
        assert!(k.contains(&FlagKind::NotFullSibling), "{flags:?}");
        for other in ["f1/C1", "f1/C2", "f1/C4", "f1/C5"] {
            assert!(kinds(&flags, other).is_empty(), "{other}: {flags:?}");
        }
    }

    #[test]
    fn ambiguous_half_sibling_names_the_candidates() {
        // Paternal haplotypes: C1 and C4 share one, C2 and C3 each have another.
        // Either C2 or C3 is the half-sibling; the data cannot tell which.
        let flags = audit(&format!("{NO_FATHER}{HALF_SIB}"));
        assert_eq!(flags.len(), 1, "{flags:?}");
        assert_eq!(flags[0].kind, FlagKind::NotFullSibling);
        assert!(flags[0].sample_id.is_none());
        assert!(flags[0].message.contains("f1/C2") && flags[0].message.contains("f1/C3"));
    }

    #[test]
    fn identical_person_in_two_families() {
        let flags = audit(&format!("{HEAD}f2,F,01:01,02:01,07:02,08:01,03:01,04:01\n"));
        assert!(
            flags
                .iter()
                .any(|f| f.kind == FlagKind::PossibleDuplicate
                    && f.sample_id.as_deref() == Some("f2/F")),
            "{flags:?}"
        );
    }
}
