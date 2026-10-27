//! Fill in alleles that typing reports leave blank, following HLA reporting
//! practice.
//!
//! * **Classical loci** (A, B, C, DRB1, DQB1, ...): a locus with one allele
//!   typed and the other blank is read as homozygous, and the typed allele is
//!   copied across. Laboratories commonly report homozygous typings this way.
//! * **DRB3, DRB4, DRB5** are present or absent depending on the DRB1
//!   haplotype, so a blank there usually means "no such gene", not
//!   "homozygous". Each DRB1 allele implies one copy of the gene of its group
//!   (DR52: DRB1*03, 11, 12, 13, 14 carry DRB3; DR53: DRB1*04, 07, 09 carry
//!   DRB4; DR51: DRB1*15, 16 carry DRB5; DRB1*01, 08, 10 carry none). Copies
//!   the DRB1 typing rules out become the absent-gene placeholder `DRBx*NNNN`;
//!   expected copies are filled from the typed alleles, a single typed allele
//!   covering two expected copies is read as homozygous, the same allele
//!   written twice where one copy is expected is read as one copy, and an expected copy
//!   with no typed allele stays unknown.
//! * A locus with neither allele typed stays untyped: nothing can be inferred.
//!
//! Filling is applied to a parsed [`Dataset`] and returns counts per locus, so
//! the choices made can be reported.

use crate::dataset::{CopiedAllele, Dataset, LocusObs};
use crate::intern::AlleleId;

/// How to read a locus with one allele typed and the other blank.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HalfTyped {
    /// The person is homozygous: copy the typed allele (reporting convention).
    Homozygous,
    /// The second allele is unknown: it may be any allele (kept as half-typed).
    Unknown,
}

impl std::str::FromStr for HalfTyped {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "homozygous" => Ok(HalfTyped::Homozygous),
            "unknown" => Ok(HalfTyped::Unknown),
            other => Err(format!(
                "unknown half-typed rule {other:?} (use homozygous or unknown)"
            )),
        }
    }
}

/// Options for [`fill_dataset`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FillOptions {
    pub half_typed: HalfTyped,
    /// Apply the DRB1-linkage rules to DRB3/4/5 (needs DRB1 in the data).
    pub drb345_from_drb1: bool,
}

impl Default for FillOptions {
    fn default() -> Self {
        FillOptions {
            half_typed: HalfTyped::Homozygous,
            drb345_from_drb1: true,
        }
    }
}

/// What was filled at one locus.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FillCount {
    pub locus: String,
    /// One typed allele copied across (homozygous).
    pub homozygous: usize,
    /// DRB3/4/5 copies set to the absent-gene placeholder `NNNN` from DRB1.
    pub absent_gene: usize,
    /// DRB3/4/5 genotypes completed from DRB1 another way (one copy filled, the
    /// other absent, or an untyped expected copy left unknown).
    pub drb_completed: usize,
    /// DRB3/4/5 written as the same allele twice where DRB1 allows one copy:
    /// read as that allele plus an absent copy (a common way of reporting).
    pub duplicate_read_as_one: usize,
    /// DRB3/4/5 typings with more copies than DRB1 allows. Left as typed and
    /// only counted: a typing error, or a rare haplotype.
    pub drb1_conflicts: usize,
}

impl FillCount {
    pub fn total(&self) -> usize {
        self.homozygous + self.absent_gene + self.drb_completed + self.duplicate_read_as_one
    }
}

/// The DRB3/4/5 gene carried with a DRB1 first field, if any.
fn drb_partner(drb1_first_field: &str) -> Option<&'static str> {
    match drb1_first_field {
        "03" | "11" | "12" | "13" | "14" => Some("DRB3"),
        "04" | "07" | "09" => Some("DRB4"),
        "15" | "16" => Some("DRB5"),
        _ => None,
    }
}

/// DRB1 first field of an allele name: `DRB1*13:02` -> `13`.
fn first_field(name: &str) -> Option<&str> {
    let typing = name.split_once('*')?.1;
    let f = typing.split(':').next()?;
    (f.len() >= 2 && f.bytes().take(2).all(|b| b.is_ascii_digit())).then(|| &f[..2])
}

/// Fill blanks in place, following [`FillOptions`]. Returns what was filled,
/// one entry per locus (zeros included).
#[allow(clippy::needless_range_loop)] // the locus index is also recorded in `ds.copied`
pub fn fill_dataset(ds: &mut Dataset, opts: &FillOptions) -> Vec<FillCount> {
    let mut counts: Vec<FillCount> = ds
        .loci
        .iter()
        .map(|l| FillCount {
            locus: l.clone(),
            ..Default::default()
        })
        .collect();
    let drb1 = ds.loci.iter().position(|l| l.eq_ignore_ascii_case("DRB1"));
    let is_drb345 = |l: &str| matches!(l.to_ascii_uppercase().as_str(), "DRB3" | "DRB4" | "DRB5");
    let rules_on = opts.drb345_from_drb1 && drb1.is_some();

    // Placeholders for absent genes, one per DRB3/4/5 locus present.
    let mut nnnn: Vec<Option<AlleleId>> = vec![None; ds.loci.len()];
    if rules_on {
        for (k, locus) in ds.loci.iter().enumerate() {
            if is_drb345(locus) {
                nnnn[k] = Some(ds.alleles.intern(k as u16, &format!("{locus}*NNNN")));
            }
        }
    }

    for s in 0..ds.samples.len() {
        // 1. Classical loci first, so a half-typed DRB1 is filled before it is used.
        // DRB3/4/5 are never copied: most people carry at most one copy of each.
        for k in 0..ds.loci.len() {
            if is_drb345(&ds.loci[k]) {
                continue;
            }
            if opts.half_typed == HalfTyped::Homozygous {
                let obs = &mut ds.samples[s].loci[k];
                if let LocusObs::OneAllele(c) = obs {
                    let typed = c.clone();
                    *obs = LocusObs::Genotypes(vec![(typed.clone(), typed.clone())]);
                    ds.copied.push(CopiedAllele {
                        sample: s,
                        locus: k,
                        typed,
                    });
                    counts[k].homozygous += 1;
                }
            }
        }
        if !rules_on {
            continue;
        }

        // 2. DRB3/4/5 from DRB1: how many copies of each gene does DRB1 imply?
        let d = drb1.unwrap();
        // Every alternative DRB1 genotype (high-resolution NGS reports list
        // several that differ only in later fields) must imply the same copies.
        let expected: Option<[usize; 3]> = match &ds.samples[s].loci[d] {
            LocusObs::Genotypes(list) if !list.is_empty() => {
                let side = |ids: &[AlleleId]| -> Option<Option<&'static str>> {
                    let mut partners = ids
                        .iter()
                        .map(|&a| first_field(ds.alleles.name(a)).map(drb_partner));
                    let first = partners.next()??;
                    // All candidates must agree on the partner gene.
                    partners.all(|p| p == Some(first)).then_some(first)
                };
                let copies = |(l, r): &(Vec<AlleleId>, Vec<AlleleId>)| -> Option<[usize; 3]> {
                    let (a, b) = (side(l)?, side(r)?);
                    let mut e = [0usize; 3];
                    for g in [a, b].into_iter().flatten() {
                        e[["DRB3", "DRB4", "DRB5"]
                            .iter()
                            .position(|x| *x == g)
                            .unwrap()] += 1;
                    }
                    Some(e)
                };
                let first = copies(&list[0]);
                // DRB1 ambiguous across groups: leave DRB3/4/5 as typed.
                first.filter(|e| list[1..].iter().all(|g| copies(g).as_ref() == Some(e)))
            }
            _ => None,
        };
        let Some(expected) = expected else { continue };
        for k in 0..ds.loci.len() {
            let Some(absent) = nnnn[k] else { continue };
            let gene = ds.loci[k].to_ascii_uppercase();
            let want = expected[["DRB3", "DRB4", "DRB5"]
                .iter()
                .position(|x| *x == gene)
                .unwrap()];
            let obs = &mut ds.samples[s].loci[k];
            let new = match (&*obs, want) {
                (LocusObs::Missing, 0) => {
                    counts[k].absent_gene += 2;
                    LocusObs::Genotypes(vec![(vec![absent], vec![absent])])
                }
                (LocusObs::Missing, 1) => {
                    counts[k].absent_gene += 1;
                    LocusObs::OneAllele(vec![absent]) // one copy absent, the other untyped
                }
                (LocusObs::OneAllele(c), 0 | 1) => {
                    counts[k].drb_completed += 1;
                    if want == 0 && c != &vec![absent] {
                        counts[k].drb1_conflicts += 1;
                    }
                    LocusObs::Genotypes(vec![(c.clone(), vec![absent])])
                }
                (LocusObs::OneAllele(c), _) => {
                    counts[k].homozygous += 1;
                    ds.copied.push(CopiedAllele {
                        sample: s,
                        locus: k,
                        typed: c.clone(),
                    });
                    LocusObs::Genotypes(vec![(c.clone(), c.clone())])
                }
                (LocusObs::Genotypes(list), 0 | 1) => {
                    let present = |side: &Vec<AlleleId>| side != &vec![absent];
                    match list.as_slice() {
                        [(l, r)] if want == 1 && l == r && present(l) => {
                            counts[k].duplicate_read_as_one += 1;
                            LocusObs::Genotypes(vec![(l.clone(), vec![absent])])
                        }
                        _ => {
                            let copies = |(l, r): &(Vec<AlleleId>, Vec<AlleleId>)| {
                                usize::from(present(l)) + usize::from(present(r))
                            };
                            if list.iter().all(|g| copies(g) > want) {
                                counts[k].drb1_conflicts += 1;
                            }
                            continue;
                        }
                    }
                }
                _ => continue, // untyped with two copies expected, or two copies typed
            };
            *obs = new;
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::{ParseOptions, read_table};

    fn load(table: &str, opts: FillOptions) -> (Dataset, Vec<FillCount>) {
        let (mut ds, _) = read_table(table.as_bytes(), &ParseOptions::default()).unwrap();
        let counts = fill_dataset(&mut ds, &opts);
        (ds, counts)
    }

    fn gl(ds: &Dataset, s: usize, locus: &str) -> String {
        let k = ds.loci.iter().position(|l| l == locus).unwrap();
        ds.samples[s].loci[k].to_gl(&ds.alleles)
    }

    #[test]
    fn half_typed_classical_locus_becomes_homozygous() {
        let table = "id,A_1,A_2,B_1,B_2\ns1,01:01,,08:01,07:02\ns2,,,08:01,08:01\n";
        let (ds, counts) = load(table, FillOptions::default());
        assert_eq!(gl(&ds, 0, "A"), "A*01:01+A*01:01");
        assert!(ds.samples[1].loci[0].is_missing(), "untyped stays untyped");
        assert_eq!(counts[0].homozygous, 1);

        let (ds, _) = load(
            table,
            FillOptions {
                half_typed: HalfTyped::Unknown,
                ..FillOptions::default()
            },
        );
        assert_eq!(gl(&ds, 0, "A"), "A*01:01+", "kept as half-typed on request");
    }

    #[test]
    fn drb345_follow_drb1() {
        // s1: DRB1*03 (DR52) + DRB1*04 (DR53): one DRB3, one DRB4, no DRB5.
        // s2: DRB1*13 + DRB1*11 (both DR52): two DRB3; one typed -> homozygous.
        // s3: DRB1*01 + DRB1*08: no DRB3/4/5 at all.
        // s4: DRB1*15 + DRB1*07: DRB5 and DRB4 expected but not typed.
        let table = "id,DRB1_1,DRB1_2,DRB3_1,DRB3_2,DRB4_1,DRB4_2,DRB5_1,DRB5_2\n\
            s1,03:01,04:01,02:02,,01:03,,,\n\
            s2,13:01,11:01,02:02,,,,,\n\
            s3,01:01,08:01,,,,,,\n\
            s4,15:01,07:01,,,,,,\n";
        let (ds, counts) = load(table, FillOptions::default());
        assert_eq!(gl(&ds, 0, "DRB3"), "DRB3*02:02+DRB3*NNNN");
        assert_eq!(gl(&ds, 0, "DRB4"), "DRB4*01:03+DRB4*NNNN");
        assert_eq!(gl(&ds, 0, "DRB5"), "DRB5*NNNN+DRB5*NNNN");
        assert_eq!(gl(&ds, 1, "DRB3"), "DRB3*02:02+DRB3*02:02");
        assert_eq!(gl(&ds, 1, "DRB4"), "DRB4*NNNN+DRB4*NNNN");
        for locus in ["DRB3", "DRB4", "DRB5"] {
            assert_eq!(gl(&ds, 2, locus), format!("{locus}*NNNN+{locus}*NNNN"));
        }
        assert_eq!(
            gl(&ds, 3, "DRB5"),
            "DRB5*NNNN+",
            "one copy expected, not typed: stays unknown"
        );
        assert_eq!(gl(&ds, 3, "DRB4"), "DRB4*NNNN+");
        assert_eq!(gl(&ds, 3, "DRB3"), "DRB3*NNNN+DRB3*NNNN");
        let drb3 = &counts[1];
        assert_eq!((drb3.homozygous, drb3.drb_completed), (1, 1));
    }

    #[test]
    fn single_drb345_copy_written_twice_is_read_as_one() {
        // s1: DRB1*01 + *11 -> one DRB3 copy, written twice.
        // s2: DRB1*03 + *11 -> two DRB3 copies, a true homozygote.
        // s3: DRB1*01 + *08 -> no DRB4 expected, but typed: a conflict, left alone.
        // s4: DRB1*01 + *11 -> one DRB3 copy, two different alleles: a conflict.
        let table = "id,DRB1_1,DRB1_2,DRB3_1,DRB3_2,DRB4_1,DRB4_2\n\
            s1,01:01,11:01,02:02,02:02,,\n\
            s2,03:01,11:01,02:02,02:02,,\n\
            s3,01:01,08:01,,,01:01,01:01\n\
            s4,01:01,11:01,02:02,01:01,,\n";
        let (ds, counts) = load(table, FillOptions::default());
        assert_eq!(gl(&ds, 0, "DRB3"), "DRB3*02:02+DRB3*NNNN");
        assert_eq!(gl(&ds, 1, "DRB3"), "DRB3*02:02+DRB3*02:02");
        assert_eq!(gl(&ds, 2, "DRB4"), "DRB4*01:01+DRB4*01:01");
        assert_eq!(gl(&ds, 3, "DRB3"), "DRB3*02:02+DRB3*01:01");
        let (drb3, drb4) = (&counts[1], &counts[2]);
        assert_eq!((drb3.duplicate_read_as_one, drb3.drb1_conflicts), (1, 1));
        assert_eq!((drb4.duplicate_read_as_one, drb4.drb1_conflicts), (0, 1));
    }

    #[test]
    fn drb345_rules_use_ambiguous_drb1_when_groups_agree() {
        // NGS reports: several DRB1 genotypes differing only in later fields.
        let input = "s1 DRB1*01:02:01:01+DRB1*03:01:01:01|DRB1*01:02:01:01+DRB1*03:01:01:05 DRB3*01:01:02:01|DRB3*01:01:02:02\n\
                     s2 DRB1*03:01+DRB1*04:01|DRB1*03:01+DRB1*15:01 DRB3*01:01\n";
        let mut ds =
            crate::parse::read_glsc(input.as_bytes(), &crate::parse::ParseOptions::default())
                .unwrap();
        fill_dataset(&mut ds, &FillOptions::default());
        let gl = |s: usize| ds.samples[s].loci[1].to_gl(&ds.alleles);
        assert_eq!(gl(0), "DRB3*01:01:02:01/DRB3*01:01:02:02+DRB3*NNNN");
        // DRB1 alternatives disagree on DR51/DR53: DRB3 left as typed.
        assert_eq!(gl(1), "DRB3*01:01+");
    }

    #[test]
    fn drb345_never_copied_without_the_drb1_rules() {
        // Without DRB1 (or with the rules off) a single DRB3 allele is kept as
        // half-typed: copying it would invent a second gene copy.
        let table = "id,DRB3_1,DRB3_2\ns1,02:02,\n";
        let (ds, _) = load(table, FillOptions::default());
        assert_eq!(gl(&ds, 0, "DRB3"), "DRB3*02:02+");
        let table = "id,DRB1_1,DRB1_2,DRB3_1,DRB3_2\ns1,03:01,04:01,02:02,\n";
        let (ds, counts) = load(
            table,
            FillOptions {
                drb345_from_drb1: false,
                ..FillOptions::default()
            },
        );
        assert_eq!(gl(&ds, 0, "DRB3"), "DRB3*02:02+");
        assert_eq!(counts.iter().map(FillCount::total).sum::<usize>(), 0);
    }

    #[test]
    fn half_typed_drb1_is_filled_before_drb345_rules() {
        // DRB1*03:01 typed once -> homozygous DR52 -> two DRB3 copies expected.
        let table = "id,DRB1_1,DRB1_2,DRB3_1,DRB3_2,DRB4_1,DRB4_2\ns1,03:01,,01:01,,,\n";
        let (ds, _) = load(table, FillOptions::default());
        assert_eq!(gl(&ds, 0, "DRB1"), "DRB1*03:01+DRB1*03:01");
        assert_eq!(gl(&ds, 0, "DRB3"), "DRB3*01:01+DRB3*01:01");
        assert_eq!(gl(&ds, 0, "DRB4"), "DRB4*NNNN+DRB4*NNNN");
    }
}
