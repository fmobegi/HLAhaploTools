//! NMDP multiple-allele codes (MACs).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::Error;

/// A large read buffer: the MAC table is ~650 MB, and small reads are very
/// slow on WSL's /mnt drives.
fn big_reader(path: impl AsRef<Path>) -> Result<BufReader<File>, Error> {
    Ok(BufReader::with_capacity(16 << 20, File::open(path)?))
}

/// Lookup table from MAC code (e.g. `AB`) to its expansion (e.g. `01/02`).
///
/// Reads either NMDP's `alpha.v3.txt` or Hapl-o-Mat's `MultipleAlleleCodes.txt`.
/// Both have one code per line, a code followed by its expansion; `alpha.v3.txt`
/// may also have a `*` column and two header lines, which are skipped.
#[derive(Debug, Default, Clone)]
pub struct MacDictionary {
    codes: HashMap<String, String>,
    /// The file's `LAST UPDATED` line, if it had one. Worth copying into
    /// output files so results can be traced back to a MAC release.
    pub source_note: Option<String>,
}

impl MacDictionary {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert(&mut self, code: &str, expansion: &str) {
        self.codes.insert(code.to_owned(), expansion.to_owned());
    }

    /// Read every code. The full NMDP table is ~1.4 million codes and ~650 MB
    /// in memory; prefer [`Self::from_reader_only`] when you know which codes
    /// an input file uses.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, Error> {
        Self::read(reader, None)
    }

    /// Read only the codes in `wanted` (see [`crate::parse::collect_mac_codes`]).
    pub fn from_reader_only<R: BufRead>(
        reader: R,
        wanted: &HashSet<String>,
    ) -> Result<Self, Error> {
        Self::read(reader, Some(wanted))
    }

    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, Error> {
        Self::from_reader(big_reader(path)?)
    }

    pub fn from_path_only(path: impl AsRef<Path>, wanted: &HashSet<String>) -> Result<Self, Error> {
        Self::from_reader_only(big_reader(path)?, wanted)
    }

    fn read<R: BufRead>(mut reader: R, wanted: Option<&HashSet<String>>) -> Result<Self, Error> {
        let mut dict = Self::new();
        let mut buf = Vec::with_capacity(1 << 16);
        loop {
            buf.clear();
            if reader.read_until(b'\n', &mut buf)? == 0 {
                break;
            }
            let Ok(line) = std::str::from_utf8(&buf) else {
                continue; // MAC files are ASCII; skip anything odd
            };
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            // NMDP writes this line in quotes: "LAST UPDATED: 10/25/26"
            let unquoted = trimmed.trim_matches('"');
            if dict.source_note.is_none() && unquoted.to_ascii_uppercase().starts_with("LAST") {
                dict.source_note = Some(unquoted.to_owned());
                continue;
            }
            let mut tokens = trimmed.split_ascii_whitespace().filter(|t| *t != "*");
            let (Some(code), Some(expansion), None) = (tokens.next(), tokens.next(), tokens.next())
            else {
                continue;
            };
            if wanted.is_some_and(|w| !w.contains(code)) {
                continue;
            }
            // Skips the "CODE SUBTYPE" header: real expansions contain digits.
            if code.bytes().all(|b| b.is_ascii_uppercase())
                && expansion.bytes().any(|b| b.is_ascii_digit())
            {
                dict.insert(code, expansion);
            }
        }
        Ok(dict)
    }

    pub fn len(&self) -> usize {
        self.codes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }

    pub fn get(&self, code: &str) -> Option<&str> {
        self.codes.get(code).map(String::as_str)
    }

    /// Expand a MAC typing such as `24:DXTW` at `locus` into full allele names.
    ///
    /// Expansion items without a colon are appended to the first field
    /// (`24:AB` with `AB = 01/02` gives `A*24:01`, `A*24:02`). Items that
    /// already contain a colon replace it (`AB = 24:02/25:01` gives
    /// `A*24:02`, `A*25:01`). This matches Hapl-o-Mat's `resolveNMDPCode`.
    pub fn expand(&self, locus: &str, typing: &str) -> Result<Vec<String>, String> {
        let (family, code) = typing
            .split_once(':')
            .ok_or_else(|| format!("{locus}*{typing} is not a MAC typing"))?;
        let expansion = self
            .get(code)
            .ok_or_else(|| format!("unknown MAC code {code} in {locus}*{typing}"))?;
        let mut out: Vec<String> = expansion
            .split('/')
            .filter(|item| !item.is_empty())
            .map(|item| {
                if item.contains(':') {
                    format!("{locus}*{item}")
                } else {
                    format!("{locus}*{family}:{item}")
                }
            })
            .collect();
        out.sort();
        out.dedup();
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALPHA_SNIPPET: &str = "\"LAST UPDATED: 10/01/26\"\n\
        *\tCODE\tSUBTYPE\n\
        \n\
        \tAB\t01/02\n\
        *\tAC\t01/03\n\
        \tDXTW\t24:02/24:03/25:01\n";

    #[test]
    fn reads_alpha_v3_layout() {
        let d = MacDictionary::from_reader(ALPHA_SNIPPET.as_bytes()).unwrap();
        assert_eq!(d.len(), 3);
        assert_eq!(d.get("AC"), Some("01/03"));
        assert_eq!(d.source_note.as_deref(), Some("LAST UPDATED: 10/01/26"));
    }

    #[test]
    fn reads_only_wanted_codes() {
        let wanted: HashSet<String> = ["AB".to_string(), "ZZ".to_string()].into();
        let d = MacDictionary::from_reader_only(ALPHA_SNIPPET.as_bytes(), &wanted).unwrap();
        assert_eq!(d.len(), 1);
        assert_eq!(d.get("AB"), Some("01/02"));
        assert_eq!(d.get("AC"), None);
    }

    #[test]
    fn expands_short_and_full_items() {
        let d = MacDictionary::from_reader(ALPHA_SNIPPET.as_bytes()).unwrap();
        assert_eq!(d.expand("A", "02:AB").unwrap(), ["A*02:01", "A*02:02"]);
        assert_eq!(
            d.expand("A", "24:DXTW").unwrap(),
            ["A*24:02", "A*24:03", "A*25:01"]
        );
        assert!(d.expand("A", "24:ZZZZ").is_err());
    }
}
