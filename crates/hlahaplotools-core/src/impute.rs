//! Per-person haplotype pairs from known haplotype frequencies (GRIMM-style imputation).
//!
//! Given a frequency table (for example `hfs.dat` from an earlier estimate) and a
//! person's typing, every pair of haplotypes in the table that fits the typing gets
//! the probability
//!
//! ```text
//! P(h1, h2 | typing) ∝ f(h1) f(h2) Σ_orders Π_loci Σ_genotypes w · fits
//! ```
//!
//! the frequency of the pair (in both orders when they differ) times the weight of the
//! person's possible genotypes it fits, normalised over all pairs. Untyped and
//! half-typed loci constrain only what is typed, so the pairs also impute them. With
//! several populations, each pair takes its haplotypes from one population, weighted
//! by the population's prior; the population posteriors come out on the way.
//!
//! Haplotypes absent from the table cannot be proposed: a person whose typing no pair
//! in the table fits gets no result and a note, as in GRIMM.

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::Error;
use crate::dataset::{Dataset, Issue};
use crate::estimate::{LocusGeno, WILD, resolve_dataset};
use crate::nomenclature::Nomenclature;
use crate::resolve::Resolution;

/// One population's haplotype frequencies.
#[derive(Debug, Clone)]
pub struct Population {
    pub name: String,
    /// `~`-joined haplotype names (loci in the dataset's order) and frequencies.
    pub haplotypes: Vec<(String, f64)>,
    /// Prior weight of the population (rescaled to sum to 1 across populations).
    pub prior: f64,
}

#[derive(Debug, Clone)]
pub struct ImputeOptions {
    /// Pairs and multilocus genotypes kept per person, most likely first.
    pub top: usize,
    /// Most candidate haplotypes considered per person and population; above it,
    /// the most frequent ones are kept and the result is marked approximate.
    pub max_candidates: usize,
}

impl Default for ImputeOptions {
    fn default() -> Self {
        ImputeOptions {
            top: 10,
            max_candidates: 5000,
        }
    }
}

/// A haplotype pair for one person.
#[derive(Debug, Clone, PartialEq)]
pub struct ImputedPair {
    pub population: String,
    pub haplotypes: [String; 2],
    pub probability: f64,
}

/// Result for one person.
#[derive(Debug, Clone)]
pub struct Imputation {
    pub sample_id: String,
    /// Most likely haplotype pairs (unordered), most likely first.
    pub pairs: Vec<ImputedPair>,
    /// Most likely multilocus unphased genotypes (`A*01:01+A*02:01^B*...`), summed over pairs.
    pub genotypes: Vec<(String, f64)>,
    /// Posterior of each population (empty with one population).
    pub populations: Vec<(String, f64)>,
    /// Number of pairs that fit the typing.
    pub n_pairs: usize,
    /// Some candidates were left out (more than `max_candidates`).
    pub approximate: bool,
}

struct Table {
    name: String,
    prior: f64,
    /// Per haplotype: code at each locus.
    keys: Vec<Vec<u32>>,
    freqs: Vec<f64>,
    names: Vec<String>,
    /// Per locus: code -> haplotypes carrying it.
    by_allele: Vec<FxHashMap<u32, Vec<u32>>>,
}

/// Impute haplotype pairs for every person in `ds`. Returns one result per person that
/// could be resolved and has at least one fitting pair, plus the reasons for the others.
pub fn impute(
    ds: &Dataset,
    nom: &Nomenclature,
    resolutions: &[Resolution],
    populations: &[Population],
    opts: &ImputeOptions,
) -> Result<(Vec<Imputation>, Vec<Issue>), Error> {
    if populations.is_empty() {
        return Err(Error::Format("no haplotype frequencies given".into()));
    }
    let (resolver, samples, mut issues) = resolve_dataset(ds, nom, resolutions)?;
    let mut codes = resolver.codes;
    let n_loci = ds.loci.len();
    let prior_sum: f64 = populations.iter().map(|p| p.prior).sum();
    let mut tables = Vec::new();
    for pop in populations {
        let mut t = Table {
            name: pop.name.clone(),
            prior: pop.prior / prior_sum,
            keys: Vec::new(),
            freqs: Vec::new(),
            names: Vec::new(),
            by_allele: vec![FxHashMap::default(); n_loci],
        };
        for (name, f) in &pop.haplotypes {
            if *f <= 0.0 {
                continue;
            }
            let parts: Vec<&str> = name.split('~').collect();
            if parts.len() != n_loci {
                return Err(Error::Format(format!(
                    "haplotype {name} in {} has {} loci, the typing has {n_loci} ({})",
                    pop.name,
                    parts.len(),
                    ds.loci.join(", ")
                )));
            }
            let mut key = Vec::with_capacity(n_loci);
            for (k, part) in parts.iter().enumerate() {
                let locus = part.split('*').next().unwrap_or("");
                if locus != ds.loci[k] {
                    return Err(Error::Format(format!(
                        "haplotype {name} in {}: locus {} where {} was expected (give the loci in \
                         the order of the frequency table)",
                        pop.name,
                        k + 1,
                        ds.loci[k]
                    )));
                }
                key.push(codes.intern(k as u16, part).0);
            }
            let h = t.keys.len() as u32;
            for (k, &c) in key.iter().enumerate() {
                t.by_allele[k].entry(c).or_default().push(h);
            }
            t.keys.push(key);
            t.freqs.push(*f);
            t.names.push(name.clone());
        }
        tables.push(t);
    }
    let codes = &codes;

    let results: Vec<Result<Imputation, Issue>> = samples
        .par_iter()
        .zip(&ds.samples)
        .filter_map(|(r, sample)| r.as_ref().map(|r| (r, sample)))
        .map(|(r, sample)| {
            let opts_k = &r.per_locus;
            let mut scored: Vec<(usize, u32, u32, f64)> = Vec::new();
            let mut approximate = false;
            for (ti, t) in tables.iter().enumerate() {
                let (cand, cut) = candidates(t, opts_k, opts.max_candidates);
                approximate |= cut;
                for (i, &x) in cand.iter().enumerate() {
                    for &y in &cand[i..] {
                        let (kx, ky) = (&t.keys[x as usize], &t.keys[y as usize]);
                        let mut s = fit(opts_k, kx, ky);
                        if x != y {
                            s += fit(opts_k, ky, kx);
                        }
                        if s > 0.0 {
                            let p = t.prior * t.freqs[x as usize] * t.freqs[y as usize] * s;
                            scored.push((ti, x, y, p));
                        }
                    }
                }
            }
            let total: f64 = scored.iter().map(|s| s.3).sum();
            if scored.is_empty() || total <= 0.0 {
                return Err(Issue {
                    line: 0,
                    sample_id: Some(sample.id.clone()),
                    message: "no pair of haplotypes in the frequency table fits this typing".into(),
                });
            }
            let mut pops = vec![0.0; tables.len()];
            let mut genos: FxHashMap<String, f64> = FxHashMap::default();
            for &(ti, x, y, p) in &scored {
                pops[ti] += p / total;
                let (kx, ky) = (&tables[ti].keys[x as usize], &tables[ti].keys[y as usize]);
                let g = kx
                    .iter()
                    .zip(ky)
                    .map(|(&a, &b)| {
                        let (a, b) = (
                            codes.name(crate::AlleleId(a)),
                            codes.name(crate::AlleleId(b)),
                        );
                        if a <= b {
                            format!("{a}+{b}")
                        } else {
                            format!("{b}+{a}")
                        }
                    })
                    .collect::<Vec<_>>()
                    .join("^");
                *genos.entry(g).or_insert(0.0) += p / total;
            }
            let n_pairs = scored.len();
            scored.sort_by(|a, b| {
                b.3.total_cmp(&a.3)
                    .then((a.0, a.1, a.2).cmp(&(b.0, b.1, b.2)))
            });
            scored.truncate(opts.top);
            let pairs = scored
                .into_iter()
                .map(|(ti, x, y, p)| ImputedPair {
                    population: tables[ti].name.clone(),
                    haplotypes: [
                        tables[ti].names[x as usize].clone(),
                        tables[ti].names[y as usize].clone(),
                    ],
                    probability: p / total,
                })
                .collect();
            let mut genotypes: Vec<(String, f64)> = genos.into_iter().collect();
            genotypes.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
            genotypes.truncate(opts.top);
            let populations = if tables.len() > 1 {
                tables.iter().map(|t| t.name.clone()).zip(pops).collect()
            } else {
                Vec::new()
            };
            Ok(Imputation {
                sample_id: sample.id.clone(),
                pairs,
                genotypes,
                populations,
                n_pairs,
                approximate,
            })
        })
        .collect();
    let mut out = Vec::new();
    for r in results {
        match r {
            Ok(x) => out.push(x),
            Err(issue) => issues.push(issue),
        }
    }
    Ok((out, issues))
}

/// Haplotypes that could be one of the person's two: at every locus, an allele the
/// person may carry (any allele where a copy is untyped). Started from the most
/// selective locus. Returns them most frequent first, and whether some were cut.
fn candidates(t: &Table, opts_k: &[Vec<(LocusGeno, f64)>], max: usize) -> (Vec<u32>, bool) {
    // Per locus: the alleles a haplotype may carry, or None for any.
    let allowed: Vec<Option<Vec<u32>>> = opts_k
        .iter()
        .map(|opts| {
            if opts.iter().any(|(g, _)| g.a == WILD || g.b == WILD) {
                return None;
            }
            let mut v: Vec<u32> = opts.iter().flat_map(|(g, _)| [g.a, g.b]).collect();
            v.sort_unstable();
            v.dedup();
            Some(v)
        })
        .collect();
    let posting_size = |k: usize, v: &[u32]| -> usize {
        v.iter()
            .map(|c| t.by_allele[k].get(c).map_or(0, Vec::len))
            .sum()
    };
    let start = allowed
        .iter()
        .enumerate()
        .filter_map(|(k, a)| a.as_ref().map(|v| (k, posting_size(k, v))))
        .min_by_key(|&(_, n)| n);
    let mut out: Vec<u32> = match start {
        Some((k, _)) => allowed[k]
            .as_ref()
            .expect("start locus is constrained")
            .iter()
            .flat_map(|c| t.by_allele[k].get(c).into_iter().flatten().copied())
            .collect(),
        None => (0..t.keys.len() as u32).collect(),
    };
    out.retain(|&h| {
        let key = &t.keys[h as usize];
        allowed
            .iter()
            .zip(key)
            .all(|(a, c)| a.as_ref().is_none_or(|v| v.binary_search(c).is_ok()))
    });
    out.sort_unstable_by(|&a, &b| {
        t.freqs[b as usize]
            .total_cmp(&t.freqs[a as usize])
            .then(a.cmp(&b))
    });
    let cut = out.len() > max;
    out.truncate(max);
    (out, cut)
}

/// Weight of the person's genotypes that `x` (first) and `y` (second) fit together.
fn fit(opts_k: &[Vec<(LocusGeno, f64)>], x: &[u32], y: &[u32]) -> f64 {
    let mut s = 1.0;
    for (k, opts) in opts_k.iter().enumerate() {
        let w: f64 = opts
            .iter()
            .filter(|(g, _)| g.allows(x[k], y[k]))
            .map(|(_, w)| w)
            .sum();
        if w <= 0.0 {
            return 0.0;
        }
        s *= w;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::{ParseOptions, read_mac};

    fn nom() -> Nomenclature {
        let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        let mut rows = Vec::new();
        for a in [
            "A*01:01", "A*02:01", "A*03:01", "B*07:02", "B*08:01", "B*44:02",
        ] {
            rows.push((a.to_string(), s(&[a])));
        }
        Nomenclature::from_rows(rows, vec![], vec![], vec![])
    }

    #[test]
    fn ranks_pairs_and_imputes_untyped_locus() {
        let text = "id\tA\tA\tB\tB\n1\t01:01\t02:01\t08:01\t07:02\n2\t01:01\t02:01\t\t\n3\t03:01\t03:01\t44:02\t44:02\n";
        let ds = read_mac(text.as_bytes(), &ParseOptions::default()).unwrap();
        let pop = Population {
            name: "all".into(),
            haplotypes: vec![
                ("A*01:01~B*08:01".into(), 0.4),
                ("A*02:01~B*07:02".into(), 0.3),
                ("A*01:01~B*07:02".into(), 0.1),
                ("A*02:01~B*08:01".into(), 0.1),
                ("A*02:01~B*44:02".into(), 0.1),
            ],
            prior: 1.0,
        };
        let r = [Resolution::TwoField; 2];
        let (out, issues) = impute(&ds, &nom(), &r, &[pop], &ImputeOptions::default()).unwrap();
        // Person 1: A*01:01~B*08:01 + A*02:01~B*07:02 (0.4*0.3) against the other phase (0.1*0.1).
        let p1 = &out[0];
        assert_eq!(p1.n_pairs, 2);
        assert_eq!(
            p1.pairs[0].haplotypes,
            ["A*01:01~B*08:01".to_string(), "A*02:01~B*07:02".into()]
        );
        assert!((p1.pairs[0].probability - 0.12 / 0.13).abs() < 1e-12);
        assert_eq!(p1.genotypes.len(), 1);
        // Person 2: B untyped, imputed; the most likely B genotype is 08:01+07:02.
        let p2 = &out[1];
        assert!(
            p2.genotypes[0].0.ends_with("B*07:02+B*08:01"),
            "{:?}",
            p2.genotypes
        );
        let total: f64 = p2.genotypes.iter().map(|g| g.1).sum();
        assert!(total <= 1.0 + 1e-12);
        // Person 3: A*03:01 is in no haplotype.
        assert_eq!(out.len(), 2);
        assert_eq!(issues.len(), 1);
        assert!(issues[0].message.contains("no pair"));
    }

    #[test]
    fn population_posteriors() {
        let text = "id\tA\tA\n1\t01:01\t01:01\n";
        let ds = read_mac(text.as_bytes(), &ParseOptions::default()).unwrap();
        let pops = [
            Population {
                name: "x".into(),
                haplotypes: vec![("A*01:01".into(), 0.5), ("A*02:01".into(), 0.5)],
                prior: 1.0,
            },
            Population {
                name: "y".into(),
                haplotypes: vec![("A*01:01".into(), 0.1), ("A*02:01".into(), 0.9)],
                prior: 1.0,
            },
        ];
        let (out, _) = impute(
            &ds,
            &nom(),
            &[Resolution::TwoField],
            &pops,
            &ImputeOptions::default(),
        )
        .unwrap();
        let px = out[0].populations[0].1;
        assert!((px - 0.25 / 0.26).abs() < 1e-12, "{px}");
        assert_eq!(out[0].pairs[0].population, "x");
    }
}
