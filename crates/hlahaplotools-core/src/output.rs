//! Result files in Hapl-o-Mat's formats.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::Error;
use crate::estimate::Estimate;

/// Haplotype frequencies: `haplotype<TAB>frequency` with 14 decimals, most
/// frequent first (Hapl-o-Mat's FILENAME_HAPLOTYPEFREQUENCIES).
pub fn write_frequencies(
    path: impl AsRef<Path>,
    est: &Estimate,
    cut: f64,
    renormalize: bool,
) -> Result<(), Error> {
    let mut out = BufWriter::new(File::create(path)?);
    for (h, f) in est.frequencies(cut, renormalize) {
        writeln!(out, "{h}\t{f:.14}")?;
    }
    out.flush()?;
    Ok(())
}

/// Convergence trace: `largest change<TAB>log-likelihood` per iteration,
/// in C-style scientific notation (Hapl-o-Mat's FILENAME_EPSILON_LOGL).
pub fn write_history(path: impl AsRef<Path>, est: &Estimate) -> Result<(), Error> {
    let mut out = BufWriter::new(File::create(path)?);
    for &(eps, log_l) in &est.history {
        writeln!(out, "{}\t{}", c_scientific(eps), c_scientific(log_l))?;
    }
    out.flush()?;
    Ok(())
}

/// Format like C's `%.14e`: `1.23456789012345e-06`.
pub fn c_scientific(v: f64) -> String {
    let s = format!("{v:.14e}");
    match s.split_once('e') {
        Some((mantissa, exp)) => {
            let (sign, digits) = exp.strip_prefix('-').map_or(("+", exp), |d| ("-", d));
            format!("{mantissa}e{sign}{digits:0>2}")
        }
        None => s, // inf / NaN
    }
}

#[cfg(test)]
mod tests {
    use super::c_scientific;

    #[test]
    fn scientific_matches_c() {
        assert_eq!(c_scientific(1e-6), "1.00000000000000e-06");
        assert_eq!(c_scientific(-1234.5), "-1.23450000000000e+03");
        assert_eq!(c_scientific(0.0), "0.00000000000000e+00");
    }
}
