//! Haplotype frequency estimation by EM, with untyped and half-typed loci.
//!
//! # Model
//!
//! Each sample is resolved to weighted multi-locus genotypes ("phenotypes"),
//! as in Hapl-o-Mat. The probability of a phenotype is the sum, over ordered
//! haplotype pairs (h, h') consistent with it, of f(h)·f(h').
//!
//! Per locus, consistency means:
//! * typed `a+b`: {h, h'} = {a, b} at that locus;
//! * half-typed `a+?`: a is one of the two alleles;
//! * untyped: anything.
//!
//! Writing "a is one of the two" as [h=a] + [h'=a] − [h=a and h'=a] and
//! multiplying over loci turns the probability into a signed sum of terms
//! F(q)·F(q'), where q and q' are *patterns*: haplotypes with some loci left
//! open (`*`), and F(q) is the summed frequency of the haplotypes matching q.
//! For fully typed samples every pattern is a single haplotype, and the
//! E and M steps reduce exactly to Hapl-o-Mat's.
//!
//! # Candidate haplotypes
//!
//! Haplotypes come from the fully specified patterns in the data. A pattern
//! with open loci that matches none of them is filled in with every allele
//! seen at those loci (up to [`EmOptions::max_completions`] per pattern), so
//! that every sample can be explained.
//!
//! # Determinism
//!
//! Work is split into fixed-size chunks that do not depend on the thread
//! count, and partial sums are merged in chunk order, so results are
//! bit-identical on any number of threads.

use rustc_hash::FxHashMap;

use rayon::prelude::*;

use crate::Error;
use crate::dataset::{Dataset, Issue, LocusObs};
use crate::intern::{AlleleId, AlleleTable};
use crate::nomenclature::Nomenclature;
use crate::resolve::{Resolution, translate};

/// Marks an open locus in a pattern.
pub(crate) const WILD: u32 = u32::MAX;
/// Hapl-o-Mat's numerical zero.
const ZERO: f64 = 1e-14;
/// Phenotypes per parallel work unit. Fixed so results don't depend on threads.
const CHUNK: usize = 512;

/// How the EM starts (Hapl-o-Mat's INITIALIZATION_HAPLOTYPEFREQUENCIES).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Init {
    /// Count how often each haplotype occurs among all phase combinations.
    NumberOccurrence,
    /// `NumberOccurrence` plus up to 10% random noise (seeded).
    Perturbation,
    /// All haplotypes equally frequent.
    Equal,
    /// Random frequencies (seeded).
    Random,
}

impl std::str::FromStr for Init {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        Ok(match s {
            "numberOccurrence" => Init::NumberOccurrence,
            "perturbation" => Init::Perturbation,
            "equal" => Init::Equal,
            "random" => Init::Random,
            other => return Err(format!("unknown initialisation {other:?}")),
        })
    }
}

#[derive(Debug, Clone)]
pub struct EmOptions {
    /// Stop when no haplotype frequency changes by more than this.
    pub epsilon: f64,
    pub max_iterations: usize,
    pub init: Init,
    pub seed: u64,
    /// Drop a sample whose typing allows more than 1/this multi-locus genotypes
    /// (Hapl-o-Mat's MINIMAL_FREQUENCY_GENOTYPES).
    pub min_genotype_frequency: f64,
    /// Most haplotypes to create when filling in one open pattern.
    pub max_completions: usize,
    /// How open patterns that match no haplotype are filled in.
    pub completion: Completion,
    /// EM runs from different start values on the same model; the best fit is kept.
    /// The first uses `init` and `seed`; the others are perturbed starts with seeds
    /// `seed + 1`, `seed + 2`, ...
    pub starts: usize,
    /// Compute each sample's most likely haplotype pair (costly with many open patterns).
    pub diplotypes: bool,
}

/// How an open pattern (a haplotype with untyped loci) that matches no haplotype
/// in the data is filled in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Completion {
    /// Each open locus takes the alleles seen at it together with the pattern's
    /// allele at its most strongly linked typed locus, in the haplotypes fully
    /// specified by the data (e.g. DQA1 given DQB1). The linked locus is the one
    /// that leaves the least uncertainty (conditional entropy) about the open
    /// locus. Falls back to every allele seen when nothing links.
    Linked,
    /// Every combination of the alleles seen at the open loci (exhaustive; can
    /// create very many haplotypes when several loci are open).
    All,
}

impl std::str::FromStr for Completion {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "linked" => Ok(Completion::Linked),
            "all" => Ok(Completion::All),
            other => Err(format!("unknown completion {other:?} (use linked or all)")),
        }
    }
}

impl Default for EmOptions {
    fn default() -> Self {
        EmOptions {
            epsilon: 1e-6,
            max_iterations: 100_000,
            init: Init::NumberOccurrence,
            seed: 0,
            min_genotype_frequency: 1e-5,
            max_completions: 10_000,
            completion: Completion::Linked,
            starts: 1,
            diplotypes: true,
        }
    }
}

/// One EM start of [`estimate`].
#[derive(Debug, Clone, PartialEq)]
pub struct StartRecord {
    pub init: Init,
    pub seed: u64,
    pub log_likelihood: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// Result of [`estimate`].
#[derive(Debug, Clone)]
pub struct Estimate {
    pub loci: Vec<String>,
    pub resolutions: Vec<Resolution>,
    /// Haplotypes with nonzero frequency, most frequent first (ties by name).
    /// Names join the per-locus codes with `~`, e.g. `A*01:01g~B*08:01g`.
    pub haplotypes: Vec<(String, f64)>,
    pub iterations: usize,
    pub converged: bool,
    pub log_likelihood: f64,
    /// Per iteration: (largest frequency change, log-likelihood after the update).
    pub history: Vec<(f64, f64)>,
    /// Samples used in the estimate.
    pub n_samples: usize,
    /// Samples used that had at least one untyped or half-typed locus.
    pub n_samples_incomplete: usize,
    pub n_phenotypes: usize,
    /// Candidate haplotypes considered (before frequencies went to zero).
    pub n_candidates: usize,
    /// Of those, how many were created to fill in open patterns.
    pub n_completions: usize,
    /// Samples left out, and why.
    pub issues: Vec<Issue>,
    /// Wall-clock seconds per stage: ("resolve", s), ("model", s), ("em", s).
    pub timings: Vec<(&'static str, f64)>,
    /// Most likely haplotype pair for each sample used (see [`Estimate::best_diplotypes`]);
    /// empty when [`EmOptions::diplotypes`] is off.
    pub diplotypes: Vec<Diplotype>,
    /// Every EM start, in the order run; the estimate is the best of them.
    pub starts: Vec<StartRecord>,
}

/// The most likely pair of haplotypes for one sample, given the estimated frequencies.
#[derive(Debug, Clone, PartialEq)]
pub struct Diplotype {
    pub sample_id: String,
    /// Haplotype names, in the same `~`-joined form as [`Estimate::haplotypes`].
    pub haplotypes: [String; 2],
    /// Posterior probability of this pair among all pairs consistent with the sample.
    pub probability: f64,
}

impl Estimate {
    /// Most likely haplotype pair per sample, in input order.
    pub fn best_diplotypes(&self) -> &[Diplotype] {
        &self.diplotypes
    }

    /// Haplotypes above `cut`, optionally rescaled to sum to 1 (as Hapl-o-Mat
    /// writes them with CUT_HAPLOTYPEFREQUENCIES / RENORMALIZE_...).
    pub fn frequencies(&self, cut: f64, renormalize: bool) -> Vec<(String, f64)> {
        let kept: Vec<_> = self
            .haplotypes
            .iter()
            .filter(|(_, f)| f - cut > ZERO)
            .cloned()
            .collect();
        if !renormalize {
            return kept;
        }
        let sum: f64 = kept.iter().map(|(_, f)| f).sum();
        kept.into_iter().map(|(h, f)| (h, f / sum)).collect()
    }
}

// ------------------------------------------------------------ resolution

/// One locus of a phenotype, as IDs in the resolved-code table.
///
/// Unphased (`phased == false`): `a+b` with a <= b, `a+?` as (a, WILD), or
/// untyped as (WILD, WILD). Phased (`phased == true`): `a` is on the first
/// haplotype and `b` on the second (either may be WILD); used when a family
/// shows which haplotype each allele is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct LocusGeno {
    pub a: u32,
    pub b: u32,
    pub phased: bool,
}

impl LocusGeno {
    pub(crate) const MISSING: LocusGeno = LocusGeno {
        a: WILD,
        b: WILD,
        phased: false,
    };

    pub(crate) fn unphased(a: u32, b: u32) -> Self {
        LocusGeno {
            a: a.min(b),
            b: a.max(b),
            phased: false,
        }
    }

    pub(crate) fn phased(a: u32, b: u32) -> Self {
        LocusGeno { a, b, phased: true }
    }

    /// Is the ordered pair (x on haplotype 1, y on haplotype 2) consistent with this locus?
    fn allows(self, x: u32, y: u32) -> bool {
        let fits = |want: u32, got: u32| want == WILD || want == got;
        if self.phased {
            return fits(self.a, x) && fits(self.b, y);
        }
        match (self.a, self.b) {
            (WILD, WILD) => true,
            (a, WILD) => x == a || y == a,
            (a, b) => (x == a && y == b) || (x == b && y == a),
        }
    }
}

/// Resolve typed alleles to target codes, caching per typed allele.
pub(crate) struct Resolver<'a> {
    ds: &'a Dataset,
    nom: &'a Nomenclature,
    resolutions: &'a [Resolution],
    pub codes: AlleleTable,
    cache: FxHashMap<AlleleId, Result<Vec<(u32, f64)>, String>>,
}

impl<'a> Resolver<'a> {
    pub(crate) fn new(
        ds: &'a Dataset,
        nom: &'a Nomenclature,
        resolutions: &'a [Resolution],
    ) -> Self {
        Resolver {
            ds,
            nom,
            resolutions,
            codes: AlleleTable::new(),
            cache: FxHashMap::default(),
        }
    }
    fn allele(&mut self, a: AlleleId) -> Result<Vec<(u32, f64)>, String> {
        if let Some(hit) = self.cache.get(&a) {
            return hit.clone();
        }
        let locus = self.ds.alleles.locus(a);
        let result = translate(
            self.nom,
            self.ds.alleles.name(a),
            self.resolutions[locus as usize],
        )
        .map(|codes| {
            codes
                .into_iter()
                .map(|(c, w)| (self.codes.intern(locus, &c).0, w))
                .collect()
        });
        self.cache.insert(a, result.clone());
        result
    }

    /// Candidates on one side share its weight equally, then each is translated.
    fn side(&mut self, candidates: &[AlleleId]) -> Result<Vec<(u32, f64)>, String> {
        let share = 1.0 / candidates.len() as f64;
        let mut out = WeightedList::default();
        for &a in candidates {
            for (code, w) in self.allele(a)? {
                out.add(code, w * share);
            }
        }
        Ok(out.items)
    }

    /// All resolved genotypes at one locus with their weights (summing to 1).
    pub(crate) fn locus(&mut self, obs: &LocusObs) -> Result<Vec<(LocusGeno, f64)>, String> {
        Ok(match obs {
            LocusObs::Missing => vec![(LocusGeno::MISSING, 1.0)],
            LocusObs::OneAllele(c) => self
                .side(c)?
                .into_iter()
                .map(|(a, w)| (LocusGeno::unphased(a, WILD), w))
                .collect(),
            LocusObs::Genotypes(list) => {
                let n = list.len() as f64;
                let mut out = WeightedList::default();
                for (left, right) in list {
                    let (l, r) = (self.side(left)?, self.side(right)?);
                    for &(a, wa) in &l {
                        for &(b, wb) in &r {
                            out.add(LocusGeno::unphased(a, b), wa * wb / n);
                        }
                    }
                }
                out.items
            }
        })
    }
}

/// Small insertion-ordered map from key to summed weight.
struct WeightedList<K> {
    items: Vec<(K, f64)>,
}

impl<K> Default for WeightedList<K> {
    fn default() -> Self {
        WeightedList { items: Vec::new() }
    }
}

impl<K: PartialEq> WeightedList<K> {
    fn add(&mut self, key: K, w: f64) {
        match self.items.iter_mut().find(|(k, _)| *k == key) {
            Some((_, total)) => *total += w,
            None => self.items.push((key, w)),
        }
    }
}

// ------------------------------------------------------------ model

/// Interned patterns: one code per locus (WILD for open loci), stored flat.
struct PatternStore {
    n_loci: usize,
    keys: Vec<u32>,
    index: FxHashMap<Box<[u32]>, u32>,
}

impl PatternStore {
    fn new(n_loci: usize) -> Self {
        PatternStore {
            n_loci,
            keys: Vec::new(),
            index: FxHashMap::default(),
        }
    }

    fn len(&self) -> usize {
        self.keys.len() / self.n_loci.max(1)
    }

    fn key(&self, p: u32) -> &[u32] {
        let start = p as usize * self.n_loci;
        &self.keys[start..start + self.n_loci]
    }

    /// Return the pattern's ID, adding it if new. Also says whether it was new.
    fn intern(&mut self, key: &[u32]) -> (u32, bool) {
        if let Some(&id) = self.index.get(key) {
            return (id, false);
        }
        let id = self.len() as u32;
        self.keys.extend_from_slice(key);
        self.index.insert(key.into(), id);
        (id, true)
    }
}

/// Everything the EM iterates over.
struct Model {
    n_loci: usize,
    /// Pattern keys, flat (n_loci codes per pattern).
    keys: Vec<u32>,
    /// Weight (expected number of samples) of each phenotype.
    num: Vec<f64>,
    /// Phenotype i's terms are terms[term_off[i]..term_off[i+1]]:
    /// (pattern x, pattern y, coefficient), x <= y, both orders merged.
    term_off: Vec<usize>,
    terms: Vec<(u32, u32, f64)>,
    /// Phenotype of each term.
    term_pheno: Vec<u32>,
    /// Where each pattern occurs: occ[occ_off[p]..occ_off[p+1]] holds
    /// term index * 2 + side (side 0: pattern is x, 1: pattern is y).
    occ_off: Vec<usize>,
    occ: Vec<u32>,
    /// For a full pattern: its haplotype index. For an open one: u32::MAX.
    pattern_hap: Vec<u32>,
    /// Open patterns and the haplotypes matching each.
    open: Vec<u32>,
    open_members: Vec<Vec<u32>>,
    /// Haplotype -> its full pattern; and -> open patterns (indices into `open`) it matches.
    hap_pattern: Vec<u32>,
    hap_open: Vec<Vec<u32>>,
    n_completions: usize,
}

/// Per-locus (x, y, sign) options for the ordered-pair expansion.
fn locus_options(g: LocusGeno) -> Vec<(u32, u32, f64)> {
    if g.phased {
        return vec![(g.a, g.b, 1.0)]; // order known from the family
    }
    match (g.a, g.b) {
        (WILD, WILD) => vec![(WILD, WILD, 1.0)],
        (a, WILD) => vec![(a, WILD, 1.0), (WILD, a, 1.0), (a, a, -1.0)],
        (a, b) if a == b => vec![(a, a, 1.0)],
        (a, b) => vec![(a, b, 1.0), (b, a, 1.0)],
    }
}

/// Visit every combination of one option per position (odometer order).
fn for_each_combination<T>(options: &[Vec<T>], mut visit: impl FnMut(&[usize])) {
    if options.iter().any(Vec::is_empty) {
        return;
    }
    let mut idx = vec![0usize; options.len()];
    loop {
        visit(&idx);
        let mut pos = options.len();
        loop {
            if pos == 0 {
                return;
            }
            pos -= 1;
            idx[pos] += 1;
            if idx[pos] < options[pos].len() {
                break;
            }
            idx[pos] = 0;
        }
    }
}

/// Loci along the MHC, telomeric class I to centromeric class II (approximate order).
const MHC_ORDER: &[&str] = &[
    "F", "G", "H", "A", "J", "E", "C", "B", "MICA", "MICB", "DRB5", "DRB4", "DRB3", "DRB1", "DQA1",
    "DQB1", "DQA2", "DQB2", "DOB", "TAP2", "TAP1", "DMB", "DMA", "DOA", "DPA1", "DPB1",
];

/// Strongly linked partners, best first: the pairs whose alleles travel together.
fn known_partners(locus: &str) -> &'static [&'static str] {
    match locus {
        "DQA1" => &["DQB1", "DRB1"],
        "DQB1" => &["DQA1", "DRB1"],
        "DRB1" => &["DQB1", "DQA1", "DRB3", "DRB4", "DRB5"],
        "DRB3" | "DRB4" | "DRB5" => &["DRB1", "DQB1"],
        "DPA1" => &["DPB1"],
        "DPB1" => &["DPA1"],
        "B" => &["C", "MICA", "MICB"],
        "C" => &["B", "A"],
        "A" => &["C", "B"],
        "MICA" => &["B", "MICB", "C"],
        "MICB" => &["MICA", "B", "C"],
        _ => &[],
    }
}

/// Smallest conditional probability P(open allele | linked allele), from the
/// two-locus estimate, for an allele to be tried when filling in.
const LINK_MIN: f64 = 1e-3;
/// Linked loci consulted per open locus (the first one that is typed is used).
const LINK_PARTNERS: usize = 3;

/// For linked completion: which alleles at an open locus go with each allele of
/// its linked loci, and which locus to consult for each open locus.
struct Linkage {
    n_loci: usize,
    /// with[l * n + m]: allele at m -> alleles at l with P(l | m) >= LINK_MIN in the
    /// two-locus EM estimate for (l, m). Only computed for the loci consulted.
    with: Vec<FxHashMap<u32, Vec<u32>>>,
    /// partners[l]: the other loci in the order they are consulted.
    partners: Vec<Vec<usize>>,
}

/// Two-locus EM estimate for loci (l, m) from every phenotype, as
/// allele at m -> alleles at l with conditional probability >= LINK_MIN.
fn pair_links(
    phenotypes: &[(Vec<LocusGeno>, f64)],
    l: usize,
    m: usize,
    loci: &[String],
    max_completions: usize,
) -> FxHashMap<u32, Vec<u32>> {
    let mut pooled: FxHashMap<Vec<LocusGeno>, f64> = FxHashMap::default();
    for (genotypes, w) in phenotypes {
        let (gl, gm) = (genotypes[l], genotypes[m]);
        if gl.a == WILD && gl.b == WILD || gm.a == WILD && gm.b == WILD {
            continue; // says nothing about how l and m go together
        }
        *pooled.entry(vec![gl, gm]).or_insert(0.0) += w;
    }
    let mut out: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    if pooled.is_empty() {
        return out;
    }
    let mut pair: Vec<(Vec<LocusGeno>, f64)> = pooled.into_iter().collect();
    pair.sort_by(|a, b| a.0.cmp(&b.0)); // deterministic order
    let n: f64 = pair.iter().map(|p| p.1).sum();
    let names = [loci[l].clone(), loci[m].clone()];
    let model = build_model(&pair, 2, max_completions, Completion::All, &names);
    let opts = EmOptions {
        epsilon: 1e-6_f64.min(0.25 / n),
        ..Default::default()
    };
    let (f, ..) = model.run(&opts, n);
    let mut by_m: FxHashMap<u32, Vec<(u32, f64)>> = FxHashMap::default();
    for (h, &v) in f.iter().enumerate() {
        if v > 0.0 {
            let key = model.key(model.hap_pattern[h]);
            by_m.entry(key[1]).or_default().push((key[0], v));
        }
    }
    for (b, list) in by_m {
        let total: f64 = list.iter().map(|x| x.1).sum();
        let mut keep: Vec<u32> = list
            .into_iter()
            .filter(|&(_, v)| v / total >= LINK_MIN)
            .map(|(a, _)| a)
            .collect();
        keep.sort_unstable();
        out.insert(b, keep);
    }
    out
}

impl Linkage {
    /// Partners: the known strongly linked loci first (DQA1-DQB1, DPA1-DPB1,
    /// DRB3/4/5-DRB1, B-C, MIC-B), then the others by distance along the MHC.
    /// Loci not on the MHC map are ordered by conditional entropy H(l | m),
    /// computed from the fully specified haplotypes. Which alleles go together is
    /// then estimated by a two-locus EM for each open locus and its first partners.
    fn new(
        keys: &[&[u32]],
        loci: &[String],
        phenotypes: &[(Vec<LocusGeno>, f64)],
        open_loci: &[bool],
        max_completions: usize,
    ) -> Self {
        let n_loci = loci.len();
        let mut with: Vec<FxHashMap<u32, Vec<u32>>> = vec![FxHashMap::default(); n_loci * n_loci];
        let mut pair_counts: Vec<FxHashMap<(u32, u32), f64>> =
            vec![FxHashMap::default(); n_loci * n_loci];
        for key in keys {
            for l in 0..n_loci {
                for m in 0..n_loci {
                    if l != m {
                        *pair_counts[l * n_loci + m]
                            .entry((key[l], key[m]))
                            .or_insert(0.0) += 1.0;
                    }
                }
            }
        }
        // H(l | m) = H(l, m) - H(m), from haplotype counts.
        let entropy = |counts: &mut dyn Iterator<Item = f64>| -> f64 {
            let v: Vec<f64> = counts.collect();
            let total: f64 = v.iter().sum();
            if total <= 0.0 {
                return 0.0;
            }
            -v.iter()
                .filter(|&&c| c > 0.0)
                .map(|&c| c / total * (c / total).ln())
                .sum::<f64>()
        };
        let mut partners = vec![Vec::new(); n_loci];
        for (l, list) in partners.iter_mut().enumerate() {
            let mut scored: Vec<(f64, usize)> = (0..n_loci)
                .filter(|&m| m != l)
                .map(|m| {
                    let joint = &pair_counts[l * n_loci + m];
                    let mut marginal: FxHashMap<u32, f64> = FxHashMap::default();
                    for (&(_, b), &c) in joint {
                        *marginal.entry(b).or_insert(0.0) += c;
                    }
                    let h = entropy(&mut joint.values().copied())
                        - entropy(&mut marginal.values().copied());
                    (h, m)
                })
                .collect();
            // Ties broken by locus order, so the choice is deterministic.
            scored.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            let by_entropy: Vec<usize> = scored.into_iter().map(|(_, m)| m).collect();
            let name = loci[l].to_ascii_uppercase();
            let pos = |m: usize| {
                MHC_ORDER
                    .iter()
                    .position(|x| *x == loci[m].to_ascii_uppercase())
            };
            let mut order: Vec<usize> = known_partners(&name)
                .iter()
                .filter_map(|p| (0..n_loci).find(|&m| m != l && loci[m].eq_ignore_ascii_case(p)))
                .collect();
            if let Some(here) = MHC_ORDER.iter().position(|x| *x == name) {
                let mut near: Vec<(usize, usize)> = (0..n_loci)
                    .filter(|&m| m != l && !order.contains(&m))
                    .filter_map(|m| pos(m).map(|p| (p.abs_diff(here), m)))
                    .collect();
                near.sort();
                order.extend(near.into_iter().map(|(_, m)| m));
            }
            let rest: Vec<usize> = by_entropy
                .into_iter()
                .filter(|m| !order.contains(m))
                .collect();
            order.extend(rest);
            *list = order;
        }
        for l in (0..n_loci).filter(|&l| open_loci[l]) {
            for &m in partners[l].iter().take(LINK_PARTNERS) {
                with[l * n_loci + m] = pair_links(phenotypes, l, m, loci, max_completions);
            }
        }
        Linkage {
            n_loci,
            with,
            partners,
        }
    }

    /// Alleles to try at open locus `l` of pattern `key`: those seen with the pattern's
    /// allele at the most informative typed locus; `None` if nothing links.
    fn allowed(&self, l: usize, key: &[u32]) -> Option<&Vec<u32>> {
        self.partners[l].iter().find_map(|&m| {
            let c = key[m];
            if c == WILD {
                return None;
            }
            self.with[l * self.n_loci + m]
                .get(&c)
                .filter(|v| !v.is_empty())
        })
    }
}

fn build_model(
    phenotypes: &[(Vec<LocusGeno>, f64)],
    n_loci: usize,
    max_completions: usize,
    completion: Completion,
    loci: &[String],
) -> Model {
    let mut store = PatternStore::new(n_loci);
    let mut num = Vec::with_capacity(phenotypes.len());
    let mut term_off = vec![0usize];
    let mut terms: Vec<(u32, u32, f64)> = Vec::new();
    let (mut x, mut y) = (Vec::with_capacity(n_loci), Vec::with_capacity(n_loci));
    let mut merged: Vec<(u32, u32, f64)> = Vec::new();

    for (genotypes, weight) in phenotypes {
        let options: Vec<_> = genotypes.iter().map(|&g| locus_options(g)).collect();
        merged.clear();
        for_each_combination(&options, |idx| {
            x.clear();
            y.clear();
            let mut coef = 1.0;
            for (locus, &i) in idx.iter().enumerate() {
                let (a, b, s) = options[locus][i];
                x.push(a);
                y.push(b);
                coef *= s;
            }
            let (px, _) = store.intern(&x);
            let (py, _) = store.intern(&y);
            let key = (px.min(py), px.max(py));
            match merged.iter_mut().find(|t| (t.0, t.1) == key) {
                Some(t) => t.2 += coef,
                None => merged.push((key.0, key.1, coef)),
            }
        });
        terms.extend(merged.iter().filter(|t| t.2 != 0.0));
        term_off.push(terms.len());
        num.push(*weight);
    }

    // Haplotypes: the full patterns, in order of first appearance.
    let n_patterns = store.len();
    let mut pattern_hap = vec![u32::MAX; n_patterns];
    let mut hap_pattern = Vec::new();
    let mut open = Vec::new();
    for p in 0..n_patterns as u32 {
        if store.key(p).contains(&WILD) {
            open.push(p);
        } else {
            pattern_hap[p as usize] = hap_pattern.len() as u32;
            hap_pattern.push(p);
        }
    }

    // Fill in open patterns that match no haplotype, using every allele
    // seen at the open loci.
    let first_members = member_lists(&store, &open, &hap_pattern);
    let mut alleles_at: Vec<Vec<u32>> = vec![Vec::new(); n_loci];
    for (genotypes, _) in phenotypes {
        for (locus, g) in genotypes.iter().enumerate() {
            alleles_at[locus].extend([g.a, g.b].into_iter().filter(|&c| c != WILD));
        }
    }
    for list in &mut alleles_at {
        list.sort_unstable();
        list.dedup();
    }
    let linkage = (completion == Completion::Linked).then(|| {
        let keys: Vec<&[u32]> = hap_pattern.iter().map(|&p| store.key(p)).collect();
        // Loci open in some pattern that will need filling in.
        let mut open_loci = vec![false; n_loci];
        for (k, &p) in open.iter().enumerate() {
            if first_members[k].is_empty() {
                for (locus, &c) in store.key(p).iter().enumerate() {
                    open_loci[locus] |= c == WILD;
                }
            }
        }
        Linkage::new(&keys, loci, phenotypes, &open_loci, max_completions)
    });
    let mut n_completions = 0;
    for (k, &p) in open.iter().enumerate() {
        if !first_members[k].is_empty() {
            continue;
        }
        let key = store.key(p).to_vec();
        let choices: Vec<Vec<u32>> = key
            .iter()
            .enumerate()
            .map(|(locus, &c)| {
                if c != WILD {
                    vec![c]
                } else if let Some(allowed) =
                    linkage.as_ref().and_then(|lk| lk.allowed(locus, &key))
                {
                    allowed.clone()
                } else {
                    alleles_at[locus].clone()
                }
            })
            .collect();
        let size = choices
            .iter()
            .try_fold(1usize, |acc, c| acc.checked_mul(c.len()));
        if size.is_none_or(|s| s > max_completions) {
            continue; // too many; samples needing it may end up unexplained
        }
        for_each_combination(&choices, |idx| {
            x.clear();
            x.extend(idx.iter().enumerate().map(|(l, &i)| choices[l][i]));
            let (id, new) = store.intern(&x);
            if new {
                pattern_hap.push(hap_pattern.len() as u32);
                hap_pattern.push(id);
                n_completions += 1;
            }
        });
    }

    let open_members = member_lists(&store, &open, &hap_pattern);
    let mut hap_open = vec![Vec::new(); hap_pattern.len()];
    for (k, list) in open_members.iter().enumerate() {
        for &h in list {
            hap_open[h as usize].push(k as u32);
        }
    }

    // Inverted index: for each pattern, the terms it appears in.
    let n_patterns = store.len();
    let mut occ_off = vec![0usize; n_patterns + 1];
    for &(px, py, _) in &terms {
        occ_off[px as usize + 1] += 1;
        if px != py {
            occ_off[py as usize + 1] += 1;
        }
    }
    for i in 0..n_patterns {
        occ_off[i + 1] += occ_off[i];
    }
    let mut fill = occ_off.clone();
    let mut occ = vec![0u32; occ_off[n_patterns]];
    let mut term_pheno = vec![0u32; terms.len()];
    for i in 0..num.len() {
        for t in term_off[i]..term_off[i + 1] {
            term_pheno[t] = i as u32;
            let (px, py, _) = terms[t];
            occ[fill[px as usize]] = (t as u32) << 1;
            fill[px as usize] += 1;
            if px != py {
                occ[fill[py as usize]] = ((t as u32) << 1) | 1;
                fill[py as usize] += 1;
            }
        }
    }

    Model {
        n_loci,
        keys: std::mem::take(&mut store.keys),
        num,
        term_off,
        terms,
        term_pheno,
        occ_off,
        occ,
        pattern_hap,
        open,
        open_members,
        hap_pattern,
        hap_open,
        n_completions,
    }
}

/// For each open pattern, the haplotypes matching it.
fn member_lists(store: &PatternStore, open: &[u32], hap_pattern: &[u32]) -> Vec<Vec<u32>> {
    let open_pos: FxHashMap<&[u32], usize> = open
        .iter()
        .enumerate()
        .map(|(k, &p)| (store.key(p), k))
        .collect();
    // Distinct sets of open loci.
    let mut masks: Vec<Vec<bool>> = open
        .iter()
        .map(|&p| store.key(p).iter().map(|&c| c == WILD).collect())
        .collect();
    masks.sort();
    masks.dedup();

    let mut lists = vec![Vec::new(); open.len()];
    let mut projected = Vec::new();
    for (h, &p) in hap_pattern.iter().enumerate() {
        let key = store.key(p);
        for mask in &masks {
            projected.clear();
            projected.extend(
                key.iter()
                    .zip(mask)
                    .map(|(&c, &m)| if m { WILD } else { c }),
            );
            if let Some(&k) = open_pos.get(projected.as_slice()) {
                lists[k].push(h as u32);
            }
        }
    }
    lists
}

// ------------------------------------------------------------ EM

/// Result of one EM run: (frequencies, iterations, converged, history, final log-likelihood).
type EmRun = (Vec<f64>, usize, bool, Vec<(f64, f64)>, f64);

/// Small deterministic RNG (SplitMix64) for the seeded initialisations.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    }
}

impl Model {
    fn n_patterns(&self) -> usize {
        self.pattern_hap.len()
    }

    fn key(&self, p: u32) -> &[u32] {
        let start = p as usize * self.n_loci;
        &self.keys[start..start + self.n_loci]
    }

    /// Probability of phenotype `i` (summed over consistent ordered haplotype pairs).
    fn phenotype_prob(&self, i: usize, big_f: &[f64]) -> f64 {
        self.terms[self.term_off[i]..self.term_off[i + 1]]
            .iter()
            .map(|&(x, y, coef)| coef * big_f[x as usize] * big_f[y as usize])
            .sum()
    }

    /// Summed frequency of every pattern.
    fn pattern_freqs(&self, f: &[f64]) -> Vec<f64> {
        let mut big_f: Vec<f64> = self
            .pattern_hap
            .par_iter()
            .map(|&h| if h == u32::MAX { 0.0 } else { f[h as usize] })
            .collect();
        let open_sums: Vec<f64> = self
            .open_members
            .par_iter()
            .map(|m| m.iter().map(|&h| f[h as usize]).sum())
            .collect();
        for (&p, s) in self.open.iter().zip(open_sums) {
            big_f[p as usize] = s;
        }
        big_f
    }

    /// One E step. Returns, per pattern, its expected count divided by its
    /// frequency (multiply by f(h) to get haplotype counts), and the
    /// log-likelihood.
    ///
    /// Each phenotype's scale num/P is computed in parallel; then each
    /// pattern gathers its own contributions in a fixed order. No step
    /// depends on how work is split between threads.
    fn expectation(&self, big_f: &[f64]) -> (Vec<f64>, f64) {
        let n_pheno = self.num.len();
        let per_chunk: Vec<(Vec<f64>, f64)> = (0..n_pheno.div_ceil(CHUNK))
            .into_par_iter()
            .map(|c| {
                let mut log_l = 0.0;
                let scale = (c * CHUNK..((c + 1) * CHUNK).min(n_pheno))
                    .map(|i| {
                        let prob: f64 = self.terms[self.term_off[i]..self.term_off[i + 1]]
                            .iter()
                            .map(|&(x, y, coef)| coef * big_f[x as usize] * big_f[y as usize])
                            .sum();
                        if prob <= 1e-300 {
                            return 0.0; // cannot be explained by the current haplotypes
                        }
                        log_l += self.num[i] * prob.ln();
                        self.num[i] / prob
                    })
                    .collect();
                (scale, log_l)
            })
            .collect();
        let mut scale = Vec::with_capacity(n_pheno);
        let mut log_l = 0.0;
        for (s, l) in per_chunk {
            scale.extend(s);
            log_l += l;
        }

        let per_pattern: Vec<f64> = (0..self.n_patterns())
            .into_par_iter()
            .map(|p| {
                let mut total = 0.0;
                for &o in &self.occ[self.occ_off[p]..self.occ_off[p + 1]] {
                    let t = (o >> 1) as usize;
                    let (x, y, coef) = self.terms[t];
                    let w = scale[self.term_pheno[t] as usize] * coef;
                    total += if x == y {
                        2.0 * w * big_f[x as usize]
                    } else if o & 1 == 0 {
                        w * big_f[y as usize]
                    } else {
                        w * big_f[x as usize]
                    };
                }
                total
            })
            .collect();
        (per_pattern, log_l)
    }

    /// Expected haplotype counts from per-pattern counts.
    fn counts(&self, f: &[f64], per_pattern: &[f64]) -> Vec<f64> {
        (0..f.len())
            .into_par_iter()
            .map(|h| {
                let mut total = per_pattern[self.hap_pattern[h] as usize];
                for &k in &self.hap_open[h] {
                    total += per_pattern[self.open[k as usize] as usize];
                }
                f[h] * total
            })
            .collect()
    }

    fn initialise(&self, init: Init, seed: u64, n_samples: f64) -> Vec<f64> {
        let n_haps = self.hap_pattern.len();
        let mut rng = SplitMix64(seed);
        match init {
            Init::Equal => return vec![1.0 / n_haps as f64; n_haps],
            Init::Random => {
                let raw: Vec<f64> = (0..n_haps).map(|_| rng.next_f64() + 1.0).collect();
                let sum: f64 = raw.iter().sum();
                return raw.into_iter().map(|v| v / sum).collect();
            }
            Init::NumberOccurrence | Init::Perturbation => {}
        }

        // Occurrence counts: each phenotype spreads its weight evenly over
        // the haplotype copies in its (positive) phase combinations.
        let mut f = vec![0.0; n_haps];
        let mut open_share = vec![0.0; self.n_patterns()];
        for i in 0..self.num.len() {
            let range = self.term_off[i]..self.term_off[i + 1];
            let total: f64 = self.terms[range.clone()]
                .iter()
                .filter(|t| t.2 > 0.0)
                .map(|t| t.2)
                .sum();
            if total <= 0.0 {
                continue;
            }
            let s = self.num[i] / (2.0 * total * n_samples);
            for &(x, y, coef) in self.terms[range].iter().filter(|t| t.2 > 0.0) {
                let pairs: &[(u32, f64)] = if x == y {
                    &[(x, 2.0)]
                } else {
                    &[(x, 1.0), (y, 1.0)]
                };
                for &(p, k) in pairs {
                    let share = k * s * coef;
                    match self.pattern_hap[p as usize] {
                        u32::MAX => open_share[p as usize] += share,
                        h => f[h as usize] += share,
                    }
                }
            }
        }
        // Shares of open patterns go to their members, in proportion to what
        // those members already have, or evenly if they have nothing yet.
        let base = f.clone();
        let mut any_open = false;
        for (k, &p) in self.open.iter().enumerate() {
            let share = open_share[p as usize];
            let members = &self.open_members[k];
            if share <= 0.0 || members.is_empty() {
                continue;
            }
            any_open = true;
            let sum: f64 = members.iter().map(|&h| base[h as usize]).sum();
            for &h in members {
                f[h as usize] += if sum > 0.0 {
                    share * base[h as usize] / sum
                } else {
                    share / members.len() as f64
                };
            }
        }
        if init == Init::Perturbation {
            for v in &mut f {
                *v *= 1.0 + rng.next_f64() / 10.0;
            }
            any_open = true;
        }
        if any_open {
            let sum: f64 = f.iter().sum();
            if sum > 0.0 {
                f.iter_mut().for_each(|v| *v /= sum);
            }
        }
        f
    }

    /// Run EM to convergence. Returns (frequencies, iterations, converged, history, final logL).
    fn run(&self, opts: &EmOptions, n_samples: f64) -> EmRun {
        let mut f = self.initialise(opts.init, opts.seed, n_samples);
        let mut history = Vec::new();
        let mut converged = false;
        let mut iterations = 0;
        let (mut per_pattern, _) = self.expectation(&self.pattern_freqs(&f));
        while iterations < opts.max_iterations {
            iterations += 1;
            let counts = self.counts(&f, &per_pattern);
            let mut largest = 0.0f64;
            for (old, c) in f.iter_mut().zip(counts) {
                let mut new = c * 0.5 / n_samples;
                if new < ZERO {
                    new = 0.0; // Hapl-o-Mat drops haplotypes this small
                }
                largest = largest.max((new - *old).abs());
                *old = new;
            }
            let (next, log_l) = self.expectation(&self.pattern_freqs(&f));
            per_pattern = next;
            history.push((largest, log_l));
            if largest - opts.epsilon < ZERO {
                converged = true;
                break;
            }
        }
        let log_l = history.last().map_or(f64::NAN, |h| h.1);
        (f, iterations, converged, history, log_l)
    }
}

// ------------------------------------------------------------ driver

/// One sample ready for estimation: per locus, its possible genotypes with weights.
#[derive(Clone)]
pub(crate) struct ResolvedSample {
    pub id: String,
    pub per_locus: Vec<Vec<(LocusGeno, f64)>>,
}

impl ResolvedSample {
    fn is_incomplete(&self) -> bool {
        self.per_locus
            .iter()
            .any(|opts| opts.iter().any(|(g, _)| g.a == WILD || g.b == WILD))
    }
}

/// The resolver (holding the code table), one entry per sample (`None` if it
/// could not be resolved), and the reasons for the `None`s.
pub(crate) type Resolved<'a> = (Resolver<'a>, Vec<Option<ResolvedSample>>, Vec<Issue>);

/// Resolve every sample of `ds` to target codes. Samples that cannot be
/// resolved (unknown alleles, nothing typed) are reported in the issue list.
pub(crate) fn resolve_dataset<'a>(
    ds: &'a Dataset,
    nom: &'a Nomenclature,
    resolutions: &'a [Resolution],
) -> Result<Resolved<'a>, Error> {
    if resolutions.len() != ds.loci.len() {
        return Err(Error::Format(format!(
            "{} resolutions given for {} loci",
            resolutions.len(),
            ds.loci.len()
        )));
    }
    let mut resolver = Resolver::new(ds, nom, resolutions);
    let mut issues = Vec::new();
    let mut out = Vec::with_capacity(ds.samples.len());
    for sample in &ds.samples {
        let issue = |message: String| Issue {
            line: 0,
            sample_id: Some(sample.id.clone()),
            message,
        };
        if sample.loci.iter().all(LocusObs::is_missing) {
            issues.push(issue("no typed loci".into()));
            out.push(None);
            continue;
        }
        match sample
            .loci
            .iter()
            .map(|obs| resolver.locus(obs))
            .collect::<Result<Vec<_>, String>>()
        {
            Ok(per_locus) => out.push(Some(ResolvedSample {
                id: sample.id.clone(),
                per_locus,
            })),
            Err(message) => {
                issues.push(issue(message));
                out.push(None);
            }
        }
    }
    Ok((resolver, out, issues))
}

/// Estimate haplotype frequencies for `ds`, resolving each locus to the given
/// resolution (one per entry of `ds.loci`).
pub fn estimate(
    ds: &Dataset,
    nom: &Nomenclature,
    resolutions: &[Resolution],
    opts: &EmOptions,
) -> Result<Estimate, Error> {
    let clock = std::time::Instant::now();
    let (resolver, samples, issues) = resolve_dataset(ds, nom, resolutions)?;
    let samples = samples.into_iter().flatten().collect();
    estimate_resolved(
        &ds.loci,
        resolutions,
        &resolver.codes,
        samples,
        issues,
        opts,
        clock,
    )
}

/// Estimate from already-resolved samples (shared by [`estimate`] and the
/// family functions).
pub(crate) fn estimate_resolved(
    loci: &[String],
    resolutions: &[Resolution],
    codes: &AlleleTable,
    samples: Vec<ResolvedSample>,
    mut issues: Vec<Issue>,
    opts: &EmOptions,
    clock: std::time::Instant,
) -> Result<Estimate, Error> {
    let n_loci = loci.len();
    let mut timings = Vec::new();
    let max_reports = 1.0 / opts.min_genotype_frequency;

    // Pool identical phenotypes; remember which ones each sample contributes to.
    let mut pheno_index: FxHashMap<Vec<LocusGeno>, usize> = FxHashMap::default();
    let mut phenotypes: Vec<(Vec<LocusGeno>, f64)> = Vec::new();
    let mut used: Vec<(String, Vec<(usize, f64)>)> = Vec::new();
    let mut n_incomplete = 0usize;
    for sample in &samples {
        // Hapl-o-Mat keeps a sample only if 1/n - min >= ZERO.
        let n_reports: f64 = sample.per_locus.iter().map(|o| o.len() as f64).product();
        if 1.0 / n_reports - opts.min_genotype_frequency < ZERO {
            issues.push(Issue {
                line: 0,
                sample_id: Some(sample.id.clone()),
                message: format!(
                    "{n_reports} possible genotypes exceed the limit of {max_reports} (min_genotype_frequency)"
                ),
            });
            continue;
        }
        if sample.is_incomplete() {
            n_incomplete += 1;
        }
        let mut reports = Vec::new();
        for_each_combination(&sample.per_locus, |idx| {
            let mut key = Vec::with_capacity(n_loci);
            let mut w = 1.0;
            for (locus, &i) in idx.iter().enumerate() {
                let (g, wi) = sample.per_locus[locus][i];
                key.push(g);
                w *= wi;
            }
            let p = match pheno_index.get(&key) {
                Some(&p) => {
                    phenotypes[p].1 += w;
                    p
                }
                None => {
                    pheno_index.insert(key.clone(), phenotypes.len());
                    phenotypes.push((key, w));
                    phenotypes.len() - 1
                }
            };
            reports.push((p, w));
        });
        used.push((sample.id.clone(), reports));
    }
    let n_samples = used.len();
    if n_samples == 0 {
        return Err(Error::Format("no usable samples".into()));
    }
    let min_epsilon = 0.5 / n_samples as f64;
    if opts.epsilon - min_epsilon > ZERO {
        return Err(Error::Format(format!(
            "epsilon {} is larger than 0.5/n = {min_epsilon}; choose a smaller epsilon",
            opts.epsilon
        )));
    }

    timings.push(("resolve", clock.elapsed().as_secs_f64()));
    let model = build_model(
        &phenotypes,
        n_loci,
        opts.max_completions,
        opts.completion,
        loci,
    );
    timings.push(("model", clock.elapsed().as_secs_f64() - timings[0].1));
    // Several starts on one model; only the best frequencies are kept.
    let mut starts = Vec::with_capacity(opts.starts.max(1));
    let mut best: Option<EmRun> = None;
    for k in 0..opts.starts.max(1) {
        let mut o = opts.clone();
        if k > 0 {
            o.init = Init::Perturbation;
            o.seed = opts.seed + k as u64;
        }
        let run = model.run(&o, n_samples as f64);
        starts.push(StartRecord {
            init: o.init,
            seed: o.seed,
            log_likelihood: run.4,
            iterations: run.1,
            converged: run.2,
        });
        // Near-ties keep the earlier start, so float noise does not decide.
        let better = match &best {
            None => true,
            Some(b) => run.4 > b.4 + 1e-9 * b.4.abs(),
        };
        if better {
            best = Some(run);
        }
    }
    let (f, iterations, converged, history, log_likelihood) = best.expect("at least one start");
    timings.push((
        "em",
        clock.elapsed().as_secs_f64() - timings[0].1 - timings[1].1,
    ));

    let hap_name = |h: u32| -> String {
        model
            .key(model.hap_pattern[h as usize])
            .iter()
            .map(|&c| codes.name(AlleleId(c)))
            .collect::<Vec<_>>()
            .join("~")
    };
    let mut haplotypes: Vec<(String, f64)> = f
        .iter()
        .enumerate()
        .filter(|(_, v)| **v > 0.0)
        .map(|(h, &v)| (hap_name(h as u32), v))
        .collect();
    haplotypes.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let diplotypes = if !opts.diplotypes {
        Vec::new()
    } else {
        best_diplotypes(&model, &f, &phenotypes, &used)
            .into_iter()
            .zip(&used)
            .map(|(best, (id, _))| match best {
                Some((h1, h2, p)) => Diplotype {
                    sample_id: id.clone(),
                    haplotypes: [hap_name(h1), hap_name(h2)],
                    probability: p,
                },
                None => Diplotype {
                    sample_id: id.clone(),
                    haplotypes: [String::new(), String::new()],
                    probability: 0.0,
                },
            })
            .collect()
    };
    timings.push((
        "diplotypes",
        clock.elapsed().as_secs_f64() - timings.iter().map(|t| t.1).sum::<f64>(),
    ));

    Ok(Estimate {
        loci: loci.to_vec(),
        resolutions: resolutions.to_vec(),
        haplotypes,
        iterations,
        converged,
        log_likelihood,
        history,
        n_samples,
        n_samples_incomplete: n_incomplete,
        n_phenotypes: phenotypes.len(),
        n_candidates: model.hap_pattern.len(),
        n_completions: model.n_completions,
        issues,
        timings,
        diplotypes,
        starts,
    })
}

/// For each sample, the most likely (unordered) haplotype pair and its
/// posterior probability.
///
/// Candidates are the best-matching haplotypes of every positive phase term;
/// each candidate's probability is then computed exactly: the frequency of
/// the pair (counting both orders when both fit) over the sample's likelihood.
fn best_diplotypes(
    model: &Model,
    f: &[f64],
    phenotypes: &[(Vec<LocusGeno>, f64)],
    used: &[(String, Vec<(usize, f64)>)],
) -> Vec<Option<(u32, u32, f64)>> {
    let big_f = model.pattern_freqs(f);
    // Best haplotype of each pattern (a full pattern is its own haplotype).
    let mut best_of = model.pattern_hap.clone();
    for (k, &p) in model.open.iter().enumerate() {
        best_of[p as usize] = model.open_members[k]
            .iter()
            .copied()
            .max_by(|&a, &b| f[a as usize].total_cmp(&f[b as usize]).then(b.cmp(&a)))
            .unwrap_or(u32::MAX);
    }
    let key_of = |h: u32| model.key(model.hap_pattern[h as usize]);
    let fits = |pheno: &[LocusGeno], x: u32, y: u32| {
        let (kx, ky) = (key_of(x), key_of(y));
        pheno
            .iter()
            .zip(kx.iter().zip(ky))
            .all(|(g, (&a, &b))| g.allows(a, b))
    };

    used.par_iter()
        .map(|(_, reports)| {
            let likelihood: f64 = reports
                .iter()
                .map(|&(p, w)| w * model.phenotype_prob(p, &big_f))
                .sum();
            if likelihood <= 0.0 {
                return None;
            }
            let mut candidates: Vec<(u32, u32)> = Vec::new();
            for &(p, _) in reports {
                for &(x, y, c) in &model.terms[model.term_off[p]..model.term_off[p + 1]] {
                    let (hx, hy) = (best_of[x as usize], best_of[y as usize]);
                    if c > 0.0 && hx != u32::MAX && hy != u32::MAX {
                        candidates.push((hx.min(hy), hx.max(hy)));
                    }
                }
            }
            candidates.sort_unstable();
            candidates.dedup();
            candidates
                .into_iter()
                .map(|(h1, h2)| {
                    let pair_f = f[h1 as usize] * f[h2 as usize];
                    let score: f64 = reports
                        .iter()
                        .map(|&(p, w)| {
                            let key = &phenotypes[p].0;
                            let orders =
                                fits(key, h1, h2) as u32 + (h1 != h2 && fits(key, h2, h1)) as u32;
                            w * pair_f * orders as f64
                        })
                        .sum();
                    (h1, h2, score / likelihood)
                })
                .max_by(|a, b| a.2.total_cmp(&b.2).then((b.0, b.1).cmp(&(a.0, a.1))))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::{ParseOptions, read_mac};

    /// Nomenclature where every typing is its own allele (asItIs makes it unused).
    fn run(input: &str, opts: &EmOptions) -> Estimate {
        let ds = read_mac(input.as_bytes(), &ParseOptions::default()).unwrap();
        let res = vec![Resolution::AsItIs; ds.loci.len()];
        estimate(&ds, &Nomenclature::default(), &res, opts).unwrap()
    }

    fn freq(e: &Estimate, h: &str) -> f64 {
        e.haplotypes.iter().find(|x| x.0 == h).map_or(0.0, |x| x.1)
    }

    #[test]
    fn homozygous_samples_give_exact_frequencies() {
        let e = run(
            "id\tA\tA\tB\tB\n1\t01:01\t01:01\t08:01\t08:01\n2\t02:01\t02:01\t07:02\t07:02\n3\t01:01\t01:01\t08:01\t08:01\n",
            &EmOptions::default(),
        );
        assert!((freq(&e, "A*01:01~B*08:01") - 2.0 / 3.0).abs() < 1e-12);
        assert!((freq(&e, "A*02:01~B*07:02") - 1.0 / 3.0).abs() < 1e-12);
        assert!(e.converged);
    }

    #[test]
    fn phase_is_resolved_by_homozygous_samples() {
        // Double heterozygote 01:01/02:01, 08:01/07:02 is phased by the
        // homozygous samples: 01:01~08:01 and 02:01~07:02 dominate.
        let mut input = String::from("id\tA\tA\tB\tB\n");
        for i in 0..20 {
            input += &format!("h{i}\t01:01\t01:01\t08:01\t08:01\n");
            input += &format!("k{i}\t02:01\t02:01\t07:02\t07:02\n");
        }
        input += "x\t01:01\t02:01\t08:01\t07:02\n";
        let e = run(
            &input,
            &EmOptions {
                epsilon: 1e-9,
                ..Default::default()
            },
        );
        assert!(freq(&e, "A*01:01~B*08:01") > 0.49);
        assert!(freq(&e, "A*01:01~B*07:02") < 1e-3);
        // The double heterozygote's most likely pair is the common phase.
        let d = e
            .best_diplotypes()
            .iter()
            .find(|d| d.sample_id == "x")
            .unwrap();
        assert_eq!(
            d.haplotypes,
            ["A*01:01~B*08:01".to_string(), "A*02:01~B*07:02".to_string()]
        );
        assert!(d.probability > 0.99, "{}", d.probability);
        // Homozygous samples are certain.
        let h = &e.best_diplotypes()[0];
        assert_eq!(h.haplotypes[0], h.haplotypes[1]);
        assert!((h.probability - 1.0).abs() < 1e-12);
    }

    #[test]
    fn untyped_locus_borrows_from_typed_samples() {
        // Samples typed at A and B fix A*01:01~B*08:01 and A*02:01~B*07:02.
        // A sample typed only at A=01:01+01:01 should not change B linkage and
        // its B alleles should follow the typed samples.
        let mut input = String::from("id\tA\tA\tB\tB\n");
        for i in 0..10 {
            input += &format!("h{i}\t01:01\t01:01\t08:01\t08:01\n");
            input += &format!("k{i}\t02:01\t02:01\t07:02\t07:02\n");
        }
        input += "m\t01:01\t01:01\t\t\n";
        let e = run(
            &input,
            &EmOptions {
                epsilon: 1e-10,
                ..Default::default()
            },
        );
        assert_eq!(e.n_samples, 21);
        assert_eq!(e.n_samples_incomplete, 1);
        assert!((freq(&e, "A*01:01~B*08:01") - 11.0 / 21.0).abs() < 1e-6);
        assert!((freq(&e, "A*02:01~B*07:02") - 10.0 / 21.0).abs() < 1e-6);
        let total: f64 = e.haplotypes.iter().map(|h| h.1).sum();
        assert!((total - 1.0).abs() < 1e-9);
    }

    #[test]
    fn half_typed_locus_counts_once() {
        // m has A=01:01 and one unknown A allele; B=08:01 homozygous.
        let mut input = String::from("id\tA\tA\tB\tB\n");
        for i in 0..10 {
            input += &format!("h{i}\t01:01\t01:01\t08:01\t08:01\n");
            input += &format!("k{i}\t02:01\t02:01\t08:01\t08:01\n");
        }
        input += "m\t01:01\t\t08:01\t08:01\n";
        let e = run(
            &input,
            &EmOptions {
                epsilon: 1e-12,
                ..Default::default()
            },
        );
        // MLE: p = freq(01:01~08:01) solves the likelihood with the
        // half-typed sample contributing p^2 + 2p(1-p) = p(2-p).
        let p = freq(&e, "A*01:01~B*08:01");
        let q = freq(&e, "A*02:01~B*08:01");
        assert!((p + q - 1.0).abs() < 1e-9);
        // Score equation: 20/p - 20/q ... solve numerically and compare.
        let score = |p: f64| 20.0 / p - 20.0 / (1.0 - p) + (2.0 - 2.0 * p) / (p * (2.0 - p));
        let (mut lo, mut hi) = (0.01, 0.99);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if score(mid) > 0.0 { lo = mid } else { hi = mid }
        }
        assert!((p - lo).abs() < 1e-6, "p = {p}, expected {lo}");
    }

    #[test]
    fn results_do_not_depend_on_thread_count() {
        let mut input = String::from("id\tA\tA\tB\tB\tC\tC\n");
        let a = ["01:01", "02:01", "03:01", "11:01", "24:02"];
        let b = ["07:02", "08:01", "44:02", "35:01", "15:01"];
        let c = ["07:01", "07:02", "04:01", "05:01", "03:04"];
        for i in 0..3000usize {
            let (x, y, z) = (i * 7 % 5, i * 3 % 5, i * 11 % 5);
            let cc = if i % 9 == 0 {
                ("", "")
            } else {
                (c[z], c[(z + i) % 5])
            };
            input += &format!(
                "s{i}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                a[x],
                a[(x + 1) % 5],
                b[y],
                b[(y + 2) % 5],
                cc.0,
                cc.1
            );
        }
        let opts = EmOptions {
            epsilon: 1e-8,
            ..Default::default()
        };
        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| run(&input, &opts));
        let many = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .unwrap()
            .install(|| run(&input, &opts));
        assert_eq!(one.iterations, many.iterations);
        assert_eq!(
            one.haplotypes, many.haplotypes,
            "bit-identical across thread counts"
        );
    }
}
