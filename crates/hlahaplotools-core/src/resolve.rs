//! Translate one typing into weighted codes at a target resolution.
//!
//! A typing like `A*02:01` at resolution g becomes a list of g groups with
//! weights summing to 1, e.g. `[("A*02:01:01G"... )]`. The rules follow
//! Hapl-o-Mat's `Allele::translateTo*` functions case by case: expand the
//! typing to full alleles, map each to the target code, and give each the
//! same share of the weight. Where the original has quirks, they are kept
//! (and marked) so results can be compared with it.

use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

use crate::nomenclature::{Lookalike, Nomenclature, cut_keeping_letter};

/// Target resolution for one locus (Hapl-o-Mat's names in brackets).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Resolution {
    /// First field, e.g. `A*01` (`1f`).
    OneField,
    /// Two fields, e.g. `A*01:01` (`2f`).
    TwoField,
    /// Three fields (`3f`).
    ThreeField,
    /// Four fields (`4f`).
    FourField,
    /// Small g groups (`g`; py-ard's `lg`).
    SmallG,
    /// Small g groups named by their two fields, without the `g` (py-ard's `lgx`):
    /// `A*01:01g` becomes `A*01:01`.
    SmallGx,
    /// G groups (`G`).
    LargeG,
    /// P groups (`P`).
    P,
    /// Leave typings exactly as written (`asItIs`).
    AsItIs,
}

impl FromStr for Resolution {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        Ok(match s {
            "1f" => Resolution::OneField,
            "2f" => Resolution::TwoField,
            "3f" => Resolution::ThreeField,
            "4f" => Resolution::FourField,
            "g" | "lg" => Resolution::SmallG,
            "lgx" => Resolution::SmallGx,
            "G" => Resolution::LargeG,
            "P" => Resolution::P,
            "asItIs" => Resolution::AsItIs,
            other => {
                return Err(format!(
                    "unknown resolution {other:?} (use 1f, 2f, 3f, 4f, g, lgx, G, P or asItIs)"
                ));
            }
        })
    }
}

impl fmt::Display for Resolution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Resolution::OneField => "1f",
            Resolution::TwoField => "2f",
            Resolution::ThreeField => "3f",
            Resolution::FourField => "4f",
            Resolution::SmallG => "g",
            Resolution::SmallGx => "lgx",
            Resolution::LargeG => "G",
            Resolution::P => "P",
            Resolution::AsItIs => "asItIs",
        })
    }
}

/// What a typing is, judged like Hapl-o-Mat's `identifyCodePrecision`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Precision {
    Fields(u8),
    SmallG,
    LargeG,
    P,
}

fn precision_of(code: &str) -> Result<Precision, String> {
    let bytes = code.as_bytes();
    let ends_with_group_letter = |letter: u8| {
        bytes.len() >= 2
            && bytes[bytes.len() - 1] == letter
            && !bytes[bytes.len() - 2].is_ascii_alphabetic()
    };
    if ends_with_group_letter(b'g') {
        return Ok(Precision::SmallG);
    }
    if ends_with_group_letter(b'G') {
        return Ok(Precision::LargeG);
    }
    if ends_with_group_letter(b'P') {
        return Ok(Precision::P);
    }
    match code.matches(':').count() {
        n @ 0..=3 => Ok(Precision::Fields(n as u8 + 1)),
        _ => Err(format!("{code} has too many fields")),
    }
}

/// Weighted codes in first-seen order; equal codes have their weights summed.
#[derive(Default)]
struct Weighted {
    codes: Vec<(String, f64)>,
    index: HashMap<String, usize>,
}

impl Weighted {
    fn add(&mut self, code: String, w: f64) {
        match self.index.get(&code) {
            Some(&i) => self.codes[i].1 += w,
            None => {
                self.index.insert(code.clone(), self.codes.len());
                self.codes.push((code, w));
            }
        }
    }

    /// Each code in `codes` (duplicates included) gets an equal share of 1.
    fn spread(codes: Vec<String>) -> Vec<(String, f64)> {
        let w = 1.0 / codes.len() as f64;
        let mut out = Weighted::default();
        for c in codes {
            out.add(c, w);
        }
        out.codes
    }
}

/// Translate one typing (full name with locus, e.g. `A*24:02`) to the target
/// resolution. Returns codes with weights summing to 1, or an error if the
/// typing is not in the nomenclature tables.
pub fn translate(
    nom: &Nomenclature,
    code: &str,
    target: Resolution,
) -> Result<Vec<(String, f64)>, String> {
    use Precision::*;
    use Resolution as R;

    if target == R::AsItIs {
        return Ok(vec![(code.to_owned(), 1.0)]);
    }
    if target == R::SmallGx {
        let mut out = Weighted::default();
        for (c, w) in translate(nom, code, R::SmallG)? {
            out.add(c.strip_suffix('g').map_or(c.clone(), str::to_owned), w);
        }
        return Ok(out.codes);
    }
    let precision = precision_of(code)?;
    if matches!(precision, Fields(_)) && nom.expand(code).is_none() {
        match nom.lookalike(code) {
            Some(Lookalike::V2(new)) if nom.read_v2 => return translate(nom, new, target),
            Some(Lookalike::V2(new)) => {
                return Err(format!(
                    "{code} is a pre-2010 allele name (now {new}); rename it, or read such \
                     names with --v2-names"
                ));
            }
            Some(Lookalike::Serology(alleles)) if nom.read_serology => {
                // Each allele carrying the antigen gets an equal share, as for a
                // first-field typing.
                let share = 1.0 / alleles.len() as f64;
                let mut out = Weighted::default();
                for a in alleles {
                    for (c, w) in translate(nom, a, target)? {
                        out.add(c, w * share);
                    }
                }
                return Ok(out.codes);
            }
            Some(Lookalike::Serology(alleles)) => {
                let (locus, antigen) = code.split_once('*').unwrap_or(("", code));
                return Err(format!(
                    "{code} is not an allele name: it looks like the serological antigen \
                     {locus}{antigen} ({} alleles). Give the allele typing, or read antigens \
                     as their alleles with --serology",
                    alleles.len()
                ));
            }
            None => {}
        }
    }

    let missing = |table: &str| {
        // A one-field typing with no alleles (`DRB1*92`) is often a serological
        // antigen entered as an allele, or a typo.
        if matches!(precision, Fields(1)) && table == "AllAllelesExpanded" {
            format!(
                "{code} is not an allele group in {table}; if it was meant as a serological \
                 antigen, give the allele typing instead (no antigen has this number either)"
            )
        } else {
            format!("{code} is not in {table} (check the IMGT/HLA release)")
        }
    };
    let expand = |c: &str| {
        nom.expand(c)
            .ok_or_else(|| missing("AllAllelesExpanded"))
            .map(<[String]>::to_vec)
    };
    let members = |c: &str| -> Result<Vec<String>, String> {
        let (table, name) = match precision {
            SmallG => (&nom.small_g, "Smallg"),
            LargeG => (&nom.large_g, "LargeG"),
            P => (&nom.p, "P"),
            Fields(_) => unreachable!("members() is only used for group typings"),
        };
        table
            .members(c)
            .map(<[String]>::to_vec)
            .ok_or_else(|| missing(name))
    };
    // Group members, each expanded to full alleles (g/P/G -> alleles).
    let members_expanded = |c: &str| -> Result<Vec<String>, String> {
        let mut out = Vec::new();
        for m in members(c)? {
            out.extend(expand(&m)?);
        }
        Ok(out)
    };
    // A typing already at the target resolution is kept as is, if known.
    let keep = |table_has: bool, table: &str| -> Result<Vec<(String, f64)>, String> {
        if table_has {
            Ok(vec![(code.to_owned(), 1.0)])
        } else {
            Err(missing(table))
        }
    };
    let cut_all = |codes: Vec<String>, n: usize| {
        codes
            .iter()
            .map(|c| cut_keeping_letter(c, n))
            .collect::<Vec<_>>()
    };

    Ok(match (target, precision) {
        // ---- n-field targets
        (R::OneField, Fields(1))
        | (R::TwoField, Fields(2))
        | (R::ThreeField, Fields(3))
        | (R::FourField, Fields(4)) => keep(nom.expand(code).is_some(), "AllAllelesExpanded")?,
        // Longer typing cut down: no expansion, weight stays whole.
        (R::OneField, Fields(_)) => {
            expand(code)?;
            vec![(cut_keeping_letter(code, 0), 1.0)]
        }
        (R::TwoField, Fields(3 | 4)) => {
            expand(code)?;
            vec![(cut_keeping_letter(code, 1), 1.0)]
        }
        (R::ThreeField, Fields(4)) => {
            expand(code)?;
            vec![(cut_keeping_letter(code, 2), 1.0)]
        }
        // Shorter typing: expand to all alleles, then cut.
        (R::TwoField, Fields(_)) => Weighted::spread(cut_all(expand(code)?, 1)),
        (R::ThreeField, Fields(_)) => Weighted::spread(cut_all(expand(code)?, 2)),
        (R::FourField, Fields(_)) => Weighted::spread(expand(code)?),
        // Groups to n fields.
        (R::OneField, _) => Weighted::spread(cut_all(members(code)?, 0)),
        (R::TwoField, _) => Weighted::spread(cut_all(members(code)?, 1)),
        (R::ThreeField, LargeG) => Weighted::spread(cut_all(members(code)?, 2)),
        (R::ThreeField, _) => Weighted::spread(cut_all(members_expanded(code)?, 2)),
        (R::FourField, _) => Weighted::spread(members_expanded(code)?),

        // ---- g target
        (R::SmallG, SmallG) => keep(nom.small_g.members(code).is_some(), "Smallg")?,
        (R::SmallG, Fields(_)) => Weighted::spread(
            expand(code)?
                .iter()
                .map(|a| nom.allele_to_small_g(a))
                .collect(),
        ),
        (R::SmallG, _) => Weighted::spread(
            members(code)?
                .iter()
                .map(|a| nom.allele_to_small_g(a))
                .collect(),
        ),

        // ---- P target
        (R::P, P) => keep(nom.p.members(code).is_some(), "P")?,
        (R::P, Fields(_)) => {
            Weighted::spread(expand(code)?.iter().map(|a| nom.allele_to_p(a)).collect())
        }
        (R::P, _) => Weighted::spread(members(code)?.iter().map(|a| nom.allele_to_p(a)).collect()),

        // ---- G target
        (R::LargeG, LargeG) => keep(nom.large_g.members(code).is_some(), "LargeG")?,
        (R::LargeG, Fields(_)) => Weighted::spread(
            expand(code)?
                .iter()
                .flat_map(|a| nom.allele_to_large_g(a))
                .collect(),
        ),
        (R::LargeG, SmallG | P) => {
            // Quirk kept from Hapl-o-Mat: each group member is looked up once
            // per allele it expands to (the inner loop never updates the code).
            let mut codes = Vec::new();
            for m in members(code)? {
                let n = expand(&m)?.len();
                for _ in 0..n {
                    codes.extend(nom.allele_to_large_g(&m));
                }
            }
            Weighted::spread(codes)
        }

        (R::AsItIs | R::SmallGx, _) => unreachable!("handled above"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    fn toy() -> Nomenclature {
        let alleles = [
            "A*01:01:01:01",
            "A*01:01:01:02N",
            "A*01:02",
            "A*02:01:01:01",
        ];
        let mut expanded: Vec<(String, Vec<String>)> = Vec::new();
        let mut add = |k: &str, v: &[&str]| expanded.push((k.to_string(), s(v)));
        add("A*01", &alleles[..3]);
        add("A*01N", &["A*01:01:01:02N"]);
        add("A*01:01", &alleles[..2]);
        add("A*01:01N", &["A*01:01:01:02N"]);
        add("A*01:02", &["A*01:02"]);
        add("A*02:01", &["A*02:01:01:01"]);
        for a in alleles {
            add(a, &[a]);
        }
        let g = vec![(
            "A*01:01:01G".to_string(),
            s(&["A*01:01:01:01", "A*01:01:01:02N"]),
        )];
        let small_g = vec![(
            "A*01:01g".to_string(),
            s(&["A*01:01:01:01", "A*01:01:01:02N", "A*01:02"]),
        )];
        let p = vec![("A*01:01P".to_string(), s(&["A*01:01:01:01", "A*01:02"]))];
        Nomenclature::from_rows(expanded, small_g, g, p)
    }

    #[test]
    fn two_field_to_g_and_p() {
        let nom = toy();
        assert_eq!(
            translate(&nom, "A*01:01", Resolution::SmallG).unwrap(),
            [("A*01:01g".into(), 1.0)]
        );
        // A*01:01:01:02N is in no P group -> cut to 2 fields keeping the N.
        assert_eq!(
            translate(&nom, "A*01:01", Resolution::P).unwrap(),
            [("A*01:01P".into(), 0.5), ("A*01:01N".into(), 0.5)]
        );
        // Not in any g group -> first two fields.
        assert_eq!(
            translate(&nom, "A*02:01", Resolution::SmallG).unwrap(),
            [("A*02:01".into(), 1.0)]
        );
    }

    #[test]
    fn first_field_expands_and_splits_weight() {
        let nom = toy();
        let got = translate(&nom, "A*01", Resolution::TwoField).unwrap();
        let total: f64 = got.iter().map(|(_, w)| w).sum();
        assert!((total - 1.0).abs() < 1e-12);
        assert_eq!(got.len(), 3); // A*01:01, A*01:01N, A*01:02
    }

    #[test]
    fn serological_antigens_and_v2_names() {
        let mut nom = toy().with_lookalikes(
            vec![
                ("A*01".into(), s(&["A*01:02"])),
                ("A*36".into(), s(&["A*01:01:01:01", "A*01:02"])),
            ],
            vec![("A*0102".into(), "A*01:02".into())],
        );
        // A known allele group is never read as an antigen.
        assert_eq!(
            translate(&nom, "A*01", Resolution::TwoField).unwrap().len(),
            3
        );
        let err = translate(&nom, "A*36", Resolution::TwoField).unwrap_err();
        assert!(
            err.contains("serological antigen A36") && err.contains("--serology"),
            "{err}"
        );
        let err = translate(&nom, "A*0102", Resolution::TwoField).unwrap_err();
        assert!(err.contains("pre-2010") && err.contains("A*01:02"), "{err}");
        let err = translate(&nom, "A*92", Resolution::TwoField).unwrap_err();
        assert!(err.contains("meant as a serological antigen"), "{err}");
        nom.read_serology = true;
        nom.read_v2 = true;
        assert_eq!(
            translate(&nom, "A*36", Resolution::TwoField).unwrap(),
            [("A*01:01".into(), 0.5), ("A*01:02".into(), 0.5)]
        );
        assert_eq!(
            translate(&nom, "A*0102", Resolution::TwoField).unwrap(),
            [("A*01:02".into(), 1.0)]
        );
    }

    #[test]
    fn groups_and_errors() {
        let nom = toy();
        assert_eq!(
            translate(&nom, "A*01:01:01G", Resolution::LargeG).unwrap(),
            [("A*01:01:01G".into(), 1.0)]
        );
        assert_eq!(
            translate(&nom, "A*01:01:01G", Resolution::SmallG).unwrap(),
            [("A*01:01g".into(), 1.0)]
        );
        assert!(translate(&nom, "A*99:01", Resolution::SmallG).is_err());
        assert_eq!(
            translate(&nom, "A*99:01", Resolution::AsItIs).unwrap(),
            [("A*99:01".into(), 1.0)]
        );
        assert_eq!("g".parse::<Resolution>(), Ok(Resolution::SmallG));
        assert_eq!("lg".parse::<Resolution>(), Ok(Resolution::SmallG));
        assert_eq!(
            translate(&nom, "A*01:01:01G", Resolution::SmallGx).unwrap(),
            [("A*01:01".into(), 1.0)]
        );
        assert!("6d".parse::<Resolution>().is_err());
    }
}
