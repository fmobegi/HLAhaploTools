//! Recognise what kind of HLA typing a string is.
//!
//! Input is the part after `LOCUS*`, e.g. `02:01:01`, `24:DXTW`, `01:01:01G`.
//! The rules follow the patterns in Hapl-o-Mat's `MAC2GLS.py`.

/// The kind of a single typing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypingKind {
    /// A named allele with 1 to 4 fields, optionally with an expression
    /// suffix (N, Q, S, A, L, C), e.g. `02:01`, `13:08Q`, `01:01:01:01`.
    Allele { fields: u8 },
    /// IMGT/HLA G group, e.g. `01:01:01G`.
    GGroup,
    /// Small g group, e.g. `01:01g`.
    SmallGGroup,
    /// P group, e.g. `01:01P`.
    PGroup,
    /// NMDP multiple-allele code, e.g. `24:DXTW` or `01:XX`.
    Mac,
}

/// Classify a typing, or return `None` if it is not valid HLA nomenclature.
pub fn classify(typing: &str) -> Option<TypingKind> {
    let last = typing.chars().last()?;

    // MAC code: one numeric field, a colon, then two or more capital letters.
    // Checked first because codes can end in G or P (`09:EMP`).
    if let Some((family, code)) = typing.split_once(':')
        && code.len() >= 2
        && code.bytes().all(|b| b.is_ascii_uppercase())
    {
        return is_field(family).then_some(TypingKind::Mac);
    }

    // Groups: numeric fields followed by G, g or P.
    if matches!(last, 'G' | 'g' | 'P') {
        let body = &typing[..typing.len() - 1];
        let n = count_numeric_fields(body)?;
        return match last {
            'G' if (2..=4).contains(&n) => Some(TypingKind::GGroup),
            'g' if n == 2 => Some(TypingKind::SmallGGroup),
            'P' if n == 2 => Some(TypingKind::PGroup),
            _ => None,
        };
    }

    // Plain allele, with an optional expression suffix.
    let body = if matches!(last, 'N' | 'Q' | 'S' | 'A' | 'L' | 'C') {
        &typing[..typing.len() - 1]
    } else {
        typing
    };
    let n = count_numeric_fields(body)?;
    (1..=4)
        .contains(&n)
        .then_some(TypingKind::Allele { fields: n as u8 })
}

/// Number of `:`-separated numeric fields, or `None` if any field is invalid.
fn count_numeric_fields(s: &str) -> Option<usize> {
    let mut n = 0;
    for field in s.split(':') {
        if !is_field(field) {
            return None;
        }
        n += 1;
    }
    Some(n)
}

/// A nomenclature field: 2 to 4 digits (`01`, `294`, `1001`).
fn is_field(s: &str) -> bool {
    (2..=4).contains(&s.len()) && s.bytes().all(|b| b.is_ascii_digit())
}

#[cfg(test)]
mod tests {
    use super::TypingKind::*;
    use super::*;

    #[test]
    fn classifies_examples_from_hapl_o_mat() {
        assert_eq!(classify("26:04"), Some(Allele { fields: 2 }));
        assert_eq!(classify("04:02:04"), Some(Allele { fields: 3 }));
        assert_eq!(classify("04:01:01:01"), Some(Allele { fields: 4 }));
        assert_eq!(classify("294:01"), Some(Allele { fields: 2 }));
        assert_eq!(classify("13:08Q"), Some(Allele { fields: 2 }));
        assert_eq!(classify("06:134N"), Some(Allele { fields: 2 }));
        assert_eq!(classify("01"), Some(Allele { fields: 1 }));
        assert_eq!(classify("24:DXTW"), Some(Mac));
        assert_eq!(classify("939:BSUFU"), Some(Mac));
        assert_eq!(classify("01:XX"), Some(Mac));
        assert_eq!(classify("09:EMP"), Some(Mac)); // ends in P, still a MAC
        assert_eq!(classify("01:ABG"), Some(Mac));
        assert_eq!(classify("15:03:01G"), Some(GGroup));
        assert_eq!(classify("26:01P"), Some(PGroup));
        assert_eq!(classify("01:01g"), Some(SmallGGroup));
    }

    #[test]
    fn rejects_garbage() {
        for bad in [
            "",
            "1",
            "01:",
            ":01",
            "01:01:01:01:01",
            "01:01X",
            "01:aB",
            "01:A",
            "AB",
            "01:01:01P",
        ] {
            assert_eq!(classify(bad), None, "{bad:?} should be rejected");
        }
    }
}
