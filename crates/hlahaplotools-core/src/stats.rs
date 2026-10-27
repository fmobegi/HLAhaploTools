//! Population-genetic statistics that go with haplotype frequencies.
//!
//! * [`hardy_weinberg`]: per locus, observed and expected heterozygosity, the
//!   inbreeding coefficient F and an exact test of Hardy–Weinberg proportions.
//!   The test compares the probability of the observed genotype table, given
//!   the allele counts, with that of tables formed by pairing the same alleles
//!   at random (Monte Carlo; seeded and independent of the thread count). Only
//!   people whose genotype at the locus is unambiguous at the target resolution
//!   are counted; untyped, half-typed and ambiguous people are counted apart.
//! * [`linkage_disequilibrium`]: per pair of loci, from two-locus haplotype
//!   frequencies estimated by EM from everyone typed at either locus: overall
//!   D′ (Hedrick 1987), Wn (Cohen 1988; Cramér's V of the haplotype table), a
//!   likelihood-ratio test of linkage equilibrium, and D′ per allele pair.
//!
//! The likelihood-ratio statistic is 2 (logL₂ − logL_A − logL_B): under linkage
//! equilibrium the two-locus genotype probability is the product of the
//! one-locus ones, so the single-locus EMs on the same people give the null
//! likelihood. It is referred to χ² with (k − 1)(l − 1) degrees of freedom,
//! which is approximate when many alleles are rare.

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::Error;
use crate::dataset::{Dataset, LocusObs, Sample};
use crate::estimate::{EmOptions, Resolver, WILD, estimate};
use crate::intern::{AlleleId, AlleleTable};
use crate::nomenclature::Nomenclature;
use crate::resolve::Resolution;

// ---------------------------------------------------------------- Hardy–Weinberg

#[derive(Debug, Clone, Copy)]
pub struct HweOptions {
    /// Random tables drawn for the exact test.
    pub permutations: usize,
    pub seed: u64,
}

impl Default for HweOptions {
    fn default() -> Self {
        HweOptions {
            permutations: 10_000,
            seed: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct HweResult {
    pub locus: String,
    /// People with an unambiguous genotype at the locus (those tested).
    pub n: usize,
    pub untyped: usize,
    pub half_typed: usize,
    /// Typing allows several genotypes at the target resolution.
    pub ambiguous: usize,
    /// Of the people tested, how many were read as homozygous from a blank second allele.
    pub homozygous_filled: usize,
    pub alleles: usize,
    pub observed_heterozygosity: f64,
    /// Nei's unbiased expected heterozygosity, 2n/(2n−1) (1 − Σ p²).
    pub expected_heterozygosity: f64,
    /// 1 − Ho / He: positive for too few heterozygotes.
    pub f: f64,
    /// Exact test: probability of a table at most as probable as the one observed.
    pub p_value: Option<f64>,
    pub permutations: usize,
}

fn ln_factorials(n: usize) -> Vec<f64> {
    let mut t = vec![0.0; n + 1];
    for i in 1..=n {
        t[i] = t[i - 1] + (i as f64).ln();
    }
    t
}

/// SplitMix64: small, seeded, reproducible.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in 0..n (n > 0), without modulo bias that matters here.
    fn below(&mut self, n: usize) -> usize {
        ((self.next() >> 11) as f64 / (1u64 << 53) as f64 * n as f64) as usize % n
    }
}

/// log P(table | allele counts) up to a constant: H ln 2 − Σ ln n_ij!.
fn table_score(
    pairs: impl Iterator<Item = (u32, u32)>,
    k: usize,
    counts: &mut [u32],
    touched: &mut Vec<usize>,
    lnf: &[f64],
) -> f64 {
    let mut het = 0usize;
    for (a, b) in pairs {
        let (i, j) = (a.min(b) as usize, a.max(b) as usize);
        het += usize::from(i != j);
        let idx = i * k + j;
        if counts[idx] == 0 {
            touched.push(idx);
        }
        counts[idx] += 1;
    }
    let mut s = het as f64 * std::f64::consts::LN_2;
    for &idx in touched.iter() {
        s -= lnf[counts[idx] as usize];
        counts[idx] = 0;
    }
    touched.clear();
    s
}

/// Exact test of Hardy–Weinberg proportions for genotypes given as allele-index pairs.
fn hwe_exact(genotypes: &[(u32, u32)], k: usize, opts: &HweOptions) -> f64 {
    let n = genotypes.len();
    let lnf = ln_factorials(n);
    let mut counts = vec![0u32; k * k];
    let mut touched = Vec::new();
    let observed = table_score(
        genotypes.iter().copied(),
        k,
        &mut counts,
        &mut touched,
        &lnf,
    );
    let tol = 1e-9 * observed.abs().max(1.0);
    let alleles: Vec<u32> = genotypes.iter().flat_map(|&(a, b)| [a, b]).collect();
    const CHUNK: usize = 500;
    let chunks = opts.permutations.div_ceil(CHUNK);
    let as_extreme: usize = (0..chunks)
        .into_par_iter()
        .map(|c| {
            // One seed per chunk: the result does not depend on the thread count.
            let mut rng = Rng(opts.seed ^ (c as u64).wrapping_mul(0xA24B_AED4_963E_E407) ^ 0x5EED);
            let mut pool = alleles.clone();
            let mut counts = vec![0u32; k * k];
            let mut touched = Vec::new();
            let todo = CHUNK.min(opts.permutations - c * CHUNK);
            let mut hits = 0usize;
            for _ in 0..todo {
                for i in (1..pool.len()).rev() {
                    pool.swap(i, rng.below(i + 1));
                }
                let s = table_score(
                    pool.as_chunks::<2>().0.iter().map(|p| (p[0], p[1])),
                    k,
                    &mut counts,
                    &mut touched,
                    &lnf,
                );
                hits += usize::from(s <= observed + tol);
            }
            hits
        })
        .sum();
    (as_extreme + 1) as f64 / (opts.permutations + 1) as f64
}

/// Hardy–Weinberg statistics for every locus of `ds` at the given resolutions.
pub fn hardy_weinberg(
    ds: &Dataset,
    nom: &Nomenclature,
    resolutions: &[Resolution],
    opts: &HweOptions,
) -> Result<Vec<HweResult>, Error> {
    if resolutions.len() != ds.loci.len() {
        return Err(Error::Format(format!(
            "{} resolutions given for {} loci",
            resolutions.len(),
            ds.loci.len()
        )));
    }
    let mut resolver = Resolver::new(ds, nom, resolutions);
    let copied: std::collections::HashSet<(usize, usize)> =
        ds.copied.iter().map(|c| (c.sample, c.locus)).collect();
    let mut out = Vec::with_capacity(ds.loci.len());
    for (l, locus) in ds.loci.iter().enumerate() {
        let (mut untyped, mut half, mut ambiguous, mut filled) = (0, 0, 0, 0);
        let mut index: FxHashMap<u32, u32> = FxHashMap::default();
        let mut genotypes: Vec<(u32, u32)> = Vec::new();
        for (s, sample) in ds.samples.iter().enumerate() {
            let obs = &sample.loci[l];
            if obs.is_missing() {
                untyped += 1;
                continue;
            }
            let Ok(options) = resolver.locus(obs) else {
                ambiguous += 1; // not resolvable at this resolution: left out
                continue;
            };
            match options.as_slice() {
                [(g, _)] if g.a != WILD && g.b != WILD => {
                    let next = index.len() as u32;
                    let a = *index.entry(g.a).or_insert(next);
                    let next = index.len() as u32;
                    let b = *index.entry(g.b).or_insert(next);
                    genotypes.push((a, b));
                    filled += usize::from(copied.contains(&(s, l)));
                }
                _ if options.iter().all(|(g, _)| g.b == WILD) => half += 1,
                _ => ambiguous += 1,
            }
        }
        let n = genotypes.len();
        let k = index.len();
        let mut allele_count = vec![0usize; k];
        for &(a, b) in &genotypes {
            allele_count[a as usize] += 1;
            allele_count[b as usize] += 1;
        }
        let het = genotypes.iter().filter(|(a, b)| a != b).count();
        let (ho, he) = if n > 0 {
            let two_n = (2 * n) as f64;
            let sum_p2: f64 = allele_count
                .iter()
                .map(|&c| (c as f64 / two_n).powi(2))
                .sum();
            let he = if n > 1 {
                two_n / (two_n - 1.0) * (1.0 - sum_p2)
            } else {
                1.0 - sum_p2
            };
            (het as f64 / n as f64, he)
        } else {
            (f64::NAN, f64::NAN)
        };
        let p_value =
            (k >= 2 && n >= 2 && opts.permutations > 0).then(|| hwe_exact(&genotypes, k, opts));
        out.push(HweResult {
            locus: locus.clone(),
            n,
            untyped,
            half_typed: half,
            ambiguous,
            homozygous_filled: filled,
            alleles: k,
            observed_heterozygosity: ho,
            expected_heterozygosity: he,
            f: if he > 0.0 { 1.0 - ho / he } else { f64::NAN },
            p_value,
            permutations: opts.permutations,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------- linkage disequilibrium

#[derive(Debug, Clone)]
pub struct LdPair {
    pub allele_a: String,
    pub allele_b: String,
    pub haplotype_frequency: f64,
    pub frequency_a: f64,
    pub frequency_b: f64,
    /// D = p_ab − p_a p_b.
    pub d: f64,
    /// D divided by its largest possible value with the same sign.
    pub d_prime: f64,
}

#[derive(Debug, Clone)]
pub struct LdResult {
    pub locus_a: String,
    pub locus_b: String,
    /// People used (typed at either locus).
    pub n: usize,
    pub alleles_a: usize,
    pub alleles_b: usize,
    /// Hedrick's overall D′ = Σ p_a p_b |D′_ab|.
    pub d_prime: f64,
    /// Cohen's Wn = sqrt(Σ D² / (p_a p_b) / (min(k, l) − 1)).
    pub wn: f64,
    /// Likelihood-ratio statistic 2 (logL₂ − logL_A − logL_B), with the nominal degrees of
    /// freedom (k − 1)(l − 1). With many rare alleles its χ² distribution does not hold, so
    /// no χ² P-value is given; use the permutation test.
    pub chi2: f64,
    pub df: usize,
    /// Permutation test of linkage equilibrium: share of datasets, with the second locus
    /// shuffled across people, whose two-locus likelihood is at least the observed one.
    /// `None` when no permutations were asked for.
    pub p_value: Option<f64>,
    pub permutations: usize,
    /// Every allele pair with a nonzero haplotype frequency, most frequent first.
    pub pairs: Vec<LdPair>,
}

/// A dataset with only the given loci (in that order), and only the given samples.
fn subset(ds: &Dataset, loci: &[usize], keep: impl Fn(&Sample) -> bool) -> Dataset {
    let mut alleles = AlleleTable::new();
    let map = |id: AlleleId, k: usize, alleles: &mut AlleleTable| {
        alleles.intern(k as u16, ds.alleles.name(id))
    };
    let mut samples = Vec::new();
    let mut kept = Vec::new();
    for (s, sample) in ds.samples.iter().enumerate() {
        if !keep(sample) {
            continue;
        }
        let per: Vec<LocusObs> = loci
            .iter()
            .enumerate()
            .map(|(k, &l)| match &sample.loci[l] {
                LocusObs::Missing => LocusObs::Missing,
                LocusObs::OneAllele(c) => {
                    LocusObs::OneAllele(c.iter().map(|&a| map(a, k, &mut alleles)).collect())
                }
                LocusObs::Genotypes(gs) => LocusObs::Genotypes(
                    gs.iter()
                        .map(|(x, y)| {
                            (
                                x.iter().map(|&a| map(a, k, &mut alleles)).collect(),
                                y.iter().map(|&a| map(a, k, &mut alleles)).collect(),
                            )
                        })
                        .collect(),
                ),
            })
            .collect();
        if per.iter().any(|o| !o.is_missing()) {
            samples.push(Sample {
                id: sample.id.clone(),
                loci: per,
            });
            kept.push(s);
        }
    }
    Dataset {
        loci: loci.iter().map(|&l| ds.loci[l].clone()).collect(),
        samples,
        alleles,
        issues: Vec::new(),
        notes: Vec::new(),
        copied: Vec::new(),
    }
}

/// Linkage disequilibrium for each pair of loci (indices into `ds.loci`).
pub fn linkage_disequilibrium(
    ds: &Dataset,
    nom: &Nomenclature,
    resolutions: &[Resolution],
    pairs: &[(usize, usize)],
    opts: &EmOptions,
    permutations: usize,
    seed: u64,
) -> Result<Vec<LdResult>, Error> {
    if resolutions.len() != ds.loci.len() {
        return Err(Error::Format(format!(
            "{} resolutions given for {} loci",
            resolutions.len(),
            ds.loci.len()
        )));
    }
    let opts = EmOptions {
        diplotypes: false,
        ..opts.clone()
    };
    let mut out = Vec::with_capacity(pairs.len());
    for &(i, j) in pairs {
        if i == j || i >= ds.loci.len() || j >= ds.loci.len() {
            return Err(Error::Format(format!("bad locus pair ({i}, {j})")));
        }
        let res2 = [resolutions[i], resolutions[j]];
        let two = subset(ds, &[i, j], |_| true);
        let est = estimate(&two, nom, &res2, &opts)?;
        // The null model must use exactly the people the two-locus estimate used.
        let left_out: std::collections::HashSet<&str> = est
            .issues
            .iter()
            .filter_map(|x| x.sample_id.as_deref())
            .collect();
        let used = |s: &Sample| !left_out.contains(s.id.as_str());
        let one_a = estimate(&subset(&two, &[0], used), nom, &res2[..1], &opts)?;
        let one_b = estimate(&subset(&two, &[1], used), nom, &res2[1..], &opts)?;
        let p_value = (permutations > 0)
            .then(|| {
                ld_permutation_p(
                    &subset(&two, &[0, 1], used),
                    nom,
                    &res2,
                    &opts,
                    est.log_likelihood,
                    permutations,
                    seed,
                )
            })
            .transpose()?;

        let mut pa: FxHashMap<&str, f64> = FxHashMap::default();
        let mut pb: FxHashMap<&str, f64> = FxHashMap::default();
        let split: Vec<(&str, &str, f64)> = est
            .haplotypes
            .iter()
            .filter_map(|(h, f)| h.split_once('~').map(|(a, b)| (a, b, *f)))
            .collect();
        for &(a, b, f) in &split {
            *pa.entry(a).or_insert(0.0) += f;
            *pb.entry(b).or_insert(0.0) += f;
        }
        let k = pa.values().filter(|&&p| p > 1e-12).count();
        let l = pb.values().filter(|&&p| p > 1e-12).count();
        let hap: FxHashMap<(&str, &str), f64> =
            split.iter().map(|&(a, b, f)| ((a, b), f)).collect();

        // Sums over every allele pair, including those with no haplotype.
        let (mut d_prime, mut chi) = (0.0, 0.0);
        let mut ld_pairs = Vec::new();
        for (&a, &p) in &pa {
            for (&b, &q) in &pb {
                if p <= 1e-12 || q <= 1e-12 {
                    continue;
                }
                let f = hap.get(&(a, b)).copied().unwrap_or(0.0);
                let d = f - p * q;
                let dmax = if d >= 0.0 {
                    (p * (1.0 - q)).min((1.0 - p) * q)
                } else {
                    (p * q).min((1.0 - p) * (1.0 - q))
                };
                let dp = if dmax > 0.0 { d / dmax } else { 0.0 };
                d_prime += p * q * dp.abs();
                chi += d * d / (p * q);
                if f > 0.0 {
                    ld_pairs.push(LdPair {
                        allele_a: a.to_owned(),
                        allele_b: b.to_owned(),
                        haplotype_frequency: f,
                        frequency_a: p,
                        frequency_b: q,
                        d,
                        d_prime: dp,
                    });
                }
            }
        }
        ld_pairs.sort_by(|x, y| {
            y.haplotype_frequency
                .total_cmp(&x.haplotype_frequency)
                .then_with(|| (&x.allele_a, &x.allele_b).cmp(&(&y.allele_a, &y.allele_b)))
        });
        let wn = if k.min(l) > 1 {
            (chi / (k.min(l) - 1) as f64).sqrt()
        } else {
            0.0
        };
        let stat =
            (2.0 * (est.log_likelihood - one_a.log_likelihood - one_b.log_likelihood)).max(0.0);
        let df = (k.saturating_sub(1)) * (l.saturating_sub(1));
        out.push(LdResult {
            locus_a: ds.loci[i].clone(),
            locus_b: ds.loci[j].clone(),
            n: est.n_samples,
            alleles_a: k,
            alleles_b: l,
            d_prime,
            wn,
            chi2: stat,
            df,
            p_value,
            permutations,
            pairs: ld_pairs,
        });
    }
    Ok(out)
}

/// Permutation test: shuffle the second locus across people (keeping each locus's own data,
/// so the single-locus likelihoods are unchanged) and re-estimate the two-locus likelihood.
fn ld_permutation_p(
    two: &Dataset,
    nom: &Nomenclature,
    res2: &[Resolution; 2],
    opts: &EmOptions,
    observed: f64,
    permutations: usize,
    seed: u64,
) -> Result<f64, Error> {
    let column: Vec<LocusObs> = two.samples.iter().map(|s| s.loci[1].clone()).collect();
    let tol = 1e-9 * observed.abs().max(1.0);
    let hits: Vec<Result<bool, Error>> = (0..permutations)
        .into_par_iter()
        .map(|r| {
            let mut rng = Rng(seed ^ (r as u64).wrapping_mul(0xD6E8_FEB8_6659_FD93) ^ 0x1D);
            let mut b = column.clone();
            for i in (1..b.len()).rev() {
                b.swap(i, rng.below(i + 1));
            }
            let mut perm = two.clone();
            for (s, obs) in perm.samples.iter_mut().zip(b) {
                s.loci[1] = obs;
            }
            perm.samples
                .retain(|s| s.loci.iter().any(|o| !o.is_missing()));
            Ok(estimate(&perm, nom, res2, opts)?.log_likelihood >= observed - tol)
        })
        .collect();
    let mut as_extreme = 0usize;
    for h in hits {
        as_extreme += usize::from(h?);
    }
    Ok((as_extreme + 1) as f64 / (permutations + 1) as f64)
}

// ---------------------------------------------------------------- χ² tail

fn ln_gamma(x: f64) -> f64 {
    // Lanczos approximation (g = 7, n = 9).
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, &c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Upper regularized incomplete gamma Q(a, x).
fn gamma_q(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    let gln = ln_gamma(a);
    if x < a + 1.0 {
        // Series for P(a, x).
        let (mut ap, mut sum, mut del) = (a, 1.0 / a, 1.0 / a);
        for _ in 0..10_000 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-15 {
                break;
            }
        }
        1.0 - (sum.ln() - x + a * x.ln() - gln).exp()
    } else {
        // Continued fraction for Q(a, x) (modified Lentz).
        let tiny = 1e-300;
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..10_000 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + an / c;
            if c.abs() < tiny {
                c = tiny;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-15 {
                break;
            }
        }
        (-x + a * x.ln() - gln).exp() * h
    }
}

/// P(χ²_df ≥ x).
pub fn chi2_sf(x: f64, df: f64) -> f64 {
    gamma_q(df / 2.0, x / 2.0).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::{ParseOptions, read_mac};

    fn data(rows: &[&str], loci: &str) -> Dataset {
        let head: Vec<String> = loci
            .split(',')
            .flat_map(|l| [l.to_owned(), l.to_owned()])
            .collect();
        let text = format!("id\t{}\n{}\n", head.join("\t"), rows.join("\n"));
        read_mac(text.as_bytes(), &ParseOptions::default()).unwrap()
    }

    #[test]
    fn chi2_tail_matches_known_values() {
        assert!((chi2_sf(3.841_458_820_694_124, 1.0) - 0.05).abs() < 1e-9);
        assert!((chi2_sf(18.307_038_053_275_14, 10.0) - 0.05).abs() < 1e-9);
        assert!((chi2_sf(0.0, 4.0) - 1.0).abs() < 1e-12);
        assert!(chi2_sf(200.0, 3.0) < 1e-40);
    }

    #[test]
    fn hwe_detects_too_few_heterozygotes_and_accepts_equilibrium() {
        // 50 x 01:01 homozygous + 50 x 02:01 homozygous: no heterozygotes at all.
        let mut rows: Vec<String> = (0..50).map(|i| format!("h{i}\t01:01\t01:01")).collect();
        rows.extend((0..50).map(|i| format!("k{i}\t02:01\t02:01")));
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let ds = data(&rows, "A");
        let r = &hardy_weinberg(
            &ds,
            &Nomenclature::default(),
            &[Resolution::AsItIs],
            &HweOptions::default(),
        )
        .unwrap()[0];
        assert_eq!((r.n, r.alleles), (100, 2));
        assert_eq!(r.observed_heterozygosity, 0.0);
        assert!(r.f > 0.99 && r.p_value.unwrap() < 1e-3);

        // Exactly Hardy–Weinberg: 25 / 50 / 25.
        let mut rows: Vec<String> = (0..25).map(|i| format!("a{i}\t01:01\t01:01")).collect();
        rows.extend((0..50).map(|i| format!("b{i}\t01:01\t02:01")));
        rows.extend((0..25).map(|i| format!("c{i}\t02:01\t02:01")));
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let ds = data(&rows, "A");
        let r = &hardy_weinberg(
            &ds,
            &Nomenclature::default(),
            &[Resolution::AsItIs],
            &HweOptions::default(),
        )
        .unwrap()[0];
        assert!(r.f.abs() < 0.02 && r.p_value.unwrap() > 0.5);
    }

    #[test]
    fn hwe_result_does_not_depend_on_threads() {
        let rows: Vec<String> = (0..300)
            .map(|i| format!("s{i}\t0{}:01\t0{}:01", 1 + i % 3, 1 + (i / 3) % 4))
            .collect();
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let ds = data(&rows, "A");
        let run = |t| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .unwrap()
                .install(|| {
                    hardy_weinberg(
                        &ds,
                        &Nomenclature::default(),
                        &[Resolution::AsItIs],
                        &HweOptions::default(),
                    )
                    .unwrap()[0]
                        .p_value
                })
        };
        assert_eq!(run(1), run(4));
    }

    #[test]
    fn ld_complete_for_linked_loci_and_absent_for_independent_ones() {
        // A*01 always with B*08, A*02 with B*07; C independent of both.
        let mut rows = Vec::new();
        for i in 0..400 {
            let (a1, b1) = if i % 2 == 0 {
                ("01:01", "08:01")
            } else {
                ("02:01", "07:02")
            };
            let (a2, b2) = if (i / 2) % 2 == 0 {
                ("01:01", "08:01")
            } else {
                ("02:01", "07:02")
            };
            let (c1, c2) = (
                ["03:04", "07:01"][(i / 4) % 2],
                ["03:04", "07:01"][(i / 8) % 2],
            );
            rows.push(format!("s{i}\t{a1}\t{a2}\t{b1}\t{b2}\t{c1}\t{c2}"));
        }
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let ds = data(&rows, "A,B,C");
        let res = [Resolution::AsItIs; 3];
        let opts = EmOptions {
            epsilon: 1e-10,
            ..EmOptions::default()
        };
        let ld = linkage_disequilibrium(
            &ds,
            &Nomenclature::default(),
            &res,
            &[(0, 1), (0, 2)],
            &opts,
            99,
            0,
        )
        .unwrap();
        let (ab, ac) = (&ld[0], &ld[1]);
        assert!(
            ab.d_prime > 0.999 && ab.wn > 0.999 && ab.p_value == Some(0.01),
            "{ab:?}"
        );
        assert!(ac.d_prime < 0.1 && ac.p_value.unwrap() > 0.05, "{ac:?}");
        assert_eq!((ab.alleles_a, ab.alleles_b, ab.df), (2, 2, 1));
        assert_eq!(ab.pairs[0].d_prime, 1.0);
    }
}
