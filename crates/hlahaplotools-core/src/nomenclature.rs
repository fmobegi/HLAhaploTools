//! Lookup tables for resolving typings: allele prefixes and g/G/P groups.
//!
//! Reads the four tables written by `scripts/build_nomenclature.py` (the same
//! formats as Hapl-o-Mat's `data/` folder, so either can be used):
//!
//! | file | line format |
//! |---|---|
//! | `AllAllelesExpanded.txt` | `prefix<TAB>allele<TAB>allele...` |
//! | `Smallg.txt`, `LargeG.txt`, `P.txt` | `group<TAB>allele<TAB>allele...` |

use std::collections::HashMap;
use std::path::Path;

use crate::Error;

/// One group table (g, G or P), keeping file order because Hapl-o-Mat's
/// lookups return the *first* matching line.
#[derive(Debug, Default, Clone)]
pub struct GroupTable {
    names: Vec<String>,
    members: Vec<Vec<String>>,
    by_name: HashMap<String, usize>,
    /// allele -> every line it appears on, in file order
    by_member: HashMap<String, Vec<usize>>,
    /// group name cut at its first `delim` -> first such line
    by_stripped_name: HashMap<String, usize>,
}

impl GroupTable {
    fn from_rows(rows: Vec<(String, Vec<String>)>, delim: Option<char>) -> Self {
        let mut t = GroupTable::default();
        for (i, (name, members)) in rows.into_iter().enumerate() {
            for m in &members {
                t.by_member.entry(m.clone()).or_default().push(i);
            }
            if let Some(d) = delim {
                // Hapl-o-Mat cuts at the *first* occurrence of the letter, so
                // for DPB1*..P the key is just "D" and never matches an allele.
                let stripped = name
                    .find(d)
                    .map_or(name.as_str(), |p| &name[..p])
                    .to_owned();
                t.by_stripped_name.entry(stripped).or_insert(i);
            }
            t.by_name.entry(name.clone()).or_insert(i);
            t.names.push(name);
            t.members.push(members);
        }
        t
    }

    /// Member alleles of a group, e.g. `A*01:01:01G`.
    pub fn members(&self, group: &str) -> Option<&[String]> {
        self.by_name.get(group).map(|&i| self.members[i].as_slice())
    }

    /// The first group (in file order) whose stripped name equals `allele`
    /// or which lists `allele` as a member.
    fn first_group_of(&self, allele: &str) -> Option<&str> {
        let by_name = self.by_stripped_name.get(allele).copied();
        let by_member = self.by_member.get(allele).and_then(|v| v.first().copied());
        let line = match (by_name, by_member) {
            (Some(a), Some(b)) => a.min(b),
            (a, b) => a.or(b)?,
        };
        Some(&self.names[line])
    }

    /// Every group (in file order) that lists `allele` as a member.
    fn all_groups_of(&self, allele: &str) -> Vec<&str> {
        self.by_member
            .get(allele)
            .map(|lines| lines.iter().map(|&i| self.names[i].as_str()).collect())
            .unwrap_or_default()
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }
}

/// All tables needed to resolve typings to a target resolution.
#[derive(Debug, Default, Clone)]
pub struct Nomenclature {
    /// typing prefix (e.g. `A*01`, `A*01:01`, `A*01:01N`) -> full allele names
    expanded: HashMap<String, Vec<String>>,
    pub small_g: GroupTable,
    pub large_g: GroupTable,
    pub p: GroupTable,
    /// Free-text label for where the tables came from, e.g. "IPD-IMGT/HLA 3.66.0".
    pub source: Option<String>,
}

fn read_rows(path: &Path) -> Result<Vec<(String, Vec<String>)>, Error> {
    // One read call: much faster than line-by-line on WSL's /mnt drives.
    let text = std::fs::read_to_string(path)
        .map_err(|e| Error::Format(format!("cannot read {}: {e}", path.display())))?;
    let mut rows = Vec::new();
    for line in text.lines() {
        let mut tokens = line.split_whitespace();
        if let Some(key) = tokens.next() {
            rows.push((key.to_owned(), tokens.map(str::to_owned).collect()));
        }
    }
    Ok(rows)
}

impl Nomenclature {
    /// Load the four tables from a folder (e.g. `data/derived/3.66.0`).
    pub fn from_dir(dir: impl AsRef<Path>) -> Result<Self, Error> {
        let dir = dir.as_ref();
        let expanded = read_rows(&dir.join("AllAllelesExpanded.txt"))?
            .into_iter()
            .collect();
        Ok(Nomenclature {
            expanded,
            small_g: GroupTable::from_rows(read_rows(&dir.join("Smallg.txt"))?, Some('g')),
            large_g: GroupTable::from_rows(read_rows(&dir.join("LargeG.txt"))?, None),
            p: GroupTable::from_rows(read_rows(&dir.join("P.txt"))?, Some('P')),
            source: dir.file_name().map(|n| n.to_string_lossy().into_owned()),
        })
    }

    /// Build from in-memory rows (used by tests).
    pub fn from_rows(
        expanded: Vec<(String, Vec<String>)>,
        small_g: Vec<(String, Vec<String>)>,
        large_g: Vec<(String, Vec<String>)>,
        p: Vec<(String, Vec<String>)>,
    ) -> Self {
        Nomenclature {
            expanded: expanded.into_iter().collect(),
            small_g: GroupTable::from_rows(small_g, Some('g')),
            large_g: GroupTable::from_rows(large_g, None),
            p: GroupTable::from_rows(p, Some('P')),
            source: None,
        }
    }

    /// Full allele names covered by a typing prefix such as `A*01:01`.
    pub fn expand(&self, code: &str) -> Option<&[String]> {
        self.expanded.get(code).map(Vec::as_slice)
    }

    /// g group of a full allele name; if none, its first two fields.
    pub fn allele_to_small_g(&self, allele: &str) -> String {
        self.small_g
            .first_group_of(allele)
            .map_or_else(|| cut_keeping_letter(allele, 1), str::to_owned)
    }

    /// P group of a full allele name; if none, its first two fields.
    pub fn allele_to_p(&self, allele: &str) -> String {
        self.p
            .first_group_of(allele)
            .map_or_else(|| cut_keeping_letter(allele, 1), str::to_owned)
    }

    /// All G groups listing a full allele name; if none, its first three fields.
    pub fn allele_to_large_g(&self, allele: &str) -> Vec<String> {
        let groups = self.large_g.all_groups_of(allele);
        if groups.is_empty() {
            vec![cut_keeping_letter(allele, 2)]
        } else {
            groups.into_iter().map(str::to_owned).collect()
        }
    }
}

/// Keep `n_colons + 1` fields; re-append a trailing expression letter if one
/// was cut off. `A*01:01:01:02N` cut to 1 colon is `A*01:01N`.
pub fn cut_keeping_letter(code: &str, n_colons: usize) -> String {
    let colons = code.matches(':').count();
    if colons <= n_colons {
        return code.to_owned();
    }
    let mut out = code;
    for _ in 0..(colons - n_colons) {
        out = &out[..out.rfind(':').expect("counted colons")];
    }
    let mut out = out.to_owned();
    if let Some(last) = code.chars().last().filter(char::is_ascii_alphabetic) {
        out.push(last);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cutting_keeps_expression_letter() {
        assert_eq!(cut_keeping_letter("A*01:01:01:02N", 1), "A*01:01N");
        assert_eq!(cut_keeping_letter("A*01:01:01:01", 0), "A*01");
        assert_eq!(cut_keeping_letter("A*01:01", 2), "A*01:01");
        // Hapl-o-Mat re-appends any final letter, including group letters.
        assert_eq!(cut_keeping_letter("A*01:01:01G", 1), "A*01:01G");
    }

    #[test]
    fn group_lookup_takes_first_line_and_dp_quirk() {
        let p = vec![
            (
                "A*01:01P".to_string(),
                vec!["A*01:01:01:01".to_string(), "A*01:32".to_string()],
            ),
            (
                "DPB1*04:01P".to_string(),
                vec!["DPB1*04:01:01:01".to_string()],
            ),
        ];
        let nom = Nomenclature::from_rows(vec![], vec![], vec![], p);
        assert_eq!(nom.allele_to_p("A*01:32"), "A*01:01P");
        // Name match: stripped "A*01:01" equals the 2-field code.
        assert_eq!(nom.allele_to_p("A*01:01"), "A*01:01P");
        // DPB1: the name is cut at the P in "DPB1", so only members match.
        assert_eq!(nom.allele_to_p("DPB1*04:01"), "DPB1*04:01");
        assert_eq!(nom.allele_to_p("DPB1*04:01:01:01"), "DPB1*04:01P");
        // Not in any group: first two fields.
        assert_eq!(nom.allele_to_p("A*02:01:01:01"), "A*02:01");
    }
}
