//! Pure-Rust core of hlahaplotools.
//!
//! Everything here is independent of Python, so `cargo test` runs without an
//! interpreter. The PyO3 bindings in the root crate are a thin layer on top.
//!
//! Pipeline so far:
//! 1. [`parse`] reads MAC or GLSC files into a [`Dataset`].
//! 2. Every allele name is interned once in an [`AlleleTable`] and from then on
//!    referred to by a compact [`AlleleId`] (`u32`).
//! 3. Each sample stores, per locus, a [`LocusObs`]: missing, one allele
//!    known, or a list of possible genotypes (unordered allele pairs).

pub mod audit;
pub mod dataset;
pub mod error;
pub mod estimate;
pub mod family;
pub mod fill;
pub mod intern;
pub mod mac;
pub mod nomenclature;
pub mod output;
pub mod parse;
pub mod resolve;
pub mod typing;

pub use audit::{AuditOptions, FlagKind, PedigreeFlag, audit_families};
pub use dataset::{Dataset, Issue, LocusObs, Sample};
pub use error::Error;
pub use estimate::{Diplotype, EmOptions, Estimate, Init, estimate};
pub use family::{
    ChildInheritance, Family, FamilyComparison, FamilyPhase, FounderComparison, MendelIssue,
    Pedigree, PhaseStatus, Role, check_mendelian, compare_to_families, estimate_families,
    phase_families, read_family_table, read_family_table_from,
};
pub use fill::{FillCount, FillOptions, HalfTyped, fill_dataset};
pub use intern::{AlleleId, AlleleTable};
pub use mac::MacDictionary;
pub use nomenclature::Nomenclature;
pub use parse::{
    ParseOptions, TableIds, collect_mac_codes, collect_mac_codes_in_file, read_glsc,
    read_glsc_file, read_mac, read_mac_file, read_table, read_table_file,
};
pub use resolve::{Resolution, translate};
pub use typing::{TypingKind, classify};
