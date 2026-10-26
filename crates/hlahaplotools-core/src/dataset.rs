use crate::intern::{AlleleId, AlleleTable};

/// One genotype as typed: candidate alleles for each of the two copies.
///
/// `A*01:01+A*02:01/A*02:02` is `([01:01], [02:01, 02:02])`. Each side's
/// candidates share that side's weight equally, as in Hapl-o-Mat.
pub type TypedGenotype = (Vec<AlleleId>, Vec<AlleleId>);

/// What is known about one sample at one locus, as typed (before resolution).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LocusObs {
    /// Locus not typed at all.
    Missing,
    /// One allele typed (one of these candidates); the other copy is unknown.
    OneAllele(Vec<AlleleId>),
    /// One or more alternative genotypes (a GL-string `|` list). Each
    /// alternative carries weight 1/n.
    Genotypes(Vec<TypedGenotype>),
}

impl LocusObs {
    pub fn is_missing(&self) -> bool {
        matches!(self, LocusObs::Missing)
    }

    /// Render as a GL string, e.g. `A*01:01+A*02:01/A*02:02`. A half-typed
    /// locus renders as `A*03:01+`; a missing one as an empty string.
    pub fn to_gl(&self, alleles: &AlleleTable) -> String {
        let list = |ids: &[AlleleId]| {
            ids.iter()
                .map(|&a| alleles.name(a))
                .collect::<Vec<_>>()
                .join("/")
        };
        match self {
            LocusObs::Missing => String::new(),
            LocusObs::OneAllele(ids) => format!("{}+", list(ids)),
            LocusObs::Genotypes(gs) => gs
                .iter()
                .map(|(l, r)| format!("{}+{}", list(l), list(r)))
                .collect::<Vec<_>>()
                .join("|"),
        }
    }
}

/// One individual: an ID plus one [`LocusObs`] per locus in [`Dataset::loci`].
#[derive(Debug, Clone)]
pub struct Sample {
    pub id: String,
    pub loci: Vec<LocusObs>,
}

/// A sample that could not be used, and why. The rest of the data is still used.
#[derive(Debug, Clone)]
pub struct Issue {
    /// 1-based line number in the input file (0 if not tied to a line).
    pub line: usize,
    pub sample_id: Option<String>,
    pub message: String,
}

/// A locus where the one typed allele was copied across (read as homozygous)
/// by [`crate::fill_dataset`]. Kept so family analyses can undo the copy when
/// it contradicts Mendelian inheritance.
#[derive(Debug, Clone)]
pub struct CopiedAllele {
    pub sample: usize,
    pub locus: usize,
    /// The typed allele (its candidates), as before the copy.
    pub typed: Vec<AlleleId>,
}

/// A parsed population sample.
#[derive(Debug, Clone, Default)]
pub struct Dataset {
    /// Locus names, in the order used by every [`Sample::loci`].
    pub loci: Vec<String>,
    pub samples: Vec<Sample>,
    /// Allele names exactly as typed (after MAC expansion).
    pub alleles: AlleleTable,
    /// Lines that could not be read.
    pub issues: Vec<Issue>,
    /// Single alleles read as homozygous by [`crate::fill_dataset`].
    pub copied: Vec<CopiedAllele>,
}

impl Dataset {
    /// Number of samples with no typing at each locus.
    pub fn missing_counts(&self) -> Vec<usize> {
        self.count_per_locus(|obs| matches!(obs, LocusObs::Missing))
    }

    /// Number of samples with only one allele typed at each locus.
    pub fn one_allele_counts(&self) -> Vec<usize> {
        self.count_per_locus(|obs| matches!(obs, LocusObs::OneAllele(_)))
    }

    fn count_per_locus(&self, pred: impl Fn(&LocusObs) -> bool) -> Vec<usize> {
        let mut counts = vec![0; self.loci.len()];
        for sample in &self.samples {
            for (count, obs) in counts.iter_mut().zip(&sample.loci) {
                if pred(obs) {
                    *count += 1;
                }
            }
        }
        counts
    }
}
