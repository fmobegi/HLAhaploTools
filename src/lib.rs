//! Python bindings. All real work happens in the `hlahaplotools-core` crate.

use pyo3::prelude::*;

#[pymodule]
mod _core {
    use std::collections::{HashMap, HashSet};
    use std::path::PathBuf;

    use hlahaplotools_core as hh;
    use pyo3::exceptions::{PyIOError, PyIndexError, PyValueError};
    use pyo3::prelude::*;
    use pyo3::types::PyDict;

    fn to_py_err(e: hh::Error) -> PyErr {
        match e {
            hh::Error::Io(e) => PyIOError::new_err(e.to_string()),
            hh::Error::Format(msg) => PyValueError::new_err(msg),
        }
    }

    fn issues_to_py(issues: &[hh::Issue]) -> Vec<(usize, Option<String>, String)> {
        issues
            .iter()
            .map(|i| (i.line, i.sample_id.clone(), i.message.clone()))
            .collect()
    }

    // ------------------------------------------------------------ MAC codes

    /// NMDP multiple-allele-code table (alpha.v3.txt or MultipleAlleleCodes.txt).
    #[pyclass(name = "MacCodes", module = "hlahaplotools._core", frozen)]
    struct MacCodes {
        inner: hh::MacDictionary,
    }

    #[pymethods]
    impl MacCodes {
        /// Load a MAC table. Pass `only` (a set of codes) to keep just those;
        /// the full NMDP table needs ~650 MB of memory.
        #[new]
        #[pyo3(signature = (path, only=None))]
        fn new(path: PathBuf, only: Option<HashSet<String>>) -> PyResult<Self> {
            let inner = match only {
                Some(codes) => hh::MacDictionary::from_path_only(path, &codes),
                None => hh::MacDictionary::from_path(path),
            };
            Ok(Self {
                inner: inner.map_err(to_py_err)?,
            })
        }

        /// The file's "LAST UPDATED" line, if present.
        #[getter]
        fn source_note(&self) -> Option<String> {
            self.inner.source_note.clone()
        }

        /// Expand e.g. ("A", "24:DXTW") into full allele names.
        fn expand(&self, locus: &str, typing: &str) -> PyResult<Vec<String>> {
            self.inner
                .expand(locus, typing)
                .map_err(PyValueError::new_err)
        }

        fn __len__(&self) -> usize {
            self.inner.len()
        }

        fn __repr__(&self) -> String {
            format!("MacCodes({} codes)", self.inner.len())
        }
    }

    // ------------------------------------------------------------ Dataset

    /// A parsed population sample.
    ///
    /// `sample(i)` gives each locus as a GL string ("A*01:01+A*02:01/A*02:02"),
    /// "A*03:01+" when only one allele is typed, or None when untyped.
    #[pyclass(name = "Dataset", module = "hlahaplotools._core", frozen)]
    struct Dataset {
        inner: hh::Dataset,
        /// What reading filled in (homozygous copies, DRB3/4/5 rules), per locus.
        fills: Vec<hh::FillCount>,
    }

    fn fill_options(half_typed: &str, drb345: bool) -> PyResult<hh::FillOptions> {
        Ok(hh::FillOptions {
            half_typed: half_typed.parse().map_err(PyValueError::new_err)?,
            drb345_from_drb1: drb345,
        })
    }

    /// {locus: {"homozygous", "absent_gene", "drb_completed", "duplicate_read_as_one", "drb1_conflicts"}},
    /// loci with fills or conflicts only.
    fn fills_to_py<'py>(py: Python<'py>, fills: &[hh::FillCount]) -> PyResult<Bound<'py, PyDict>> {
        let out = PyDict::new(py);
        for f in fills.iter().filter(|f| f.total() + f.drb1_conflicts > 0) {
            let d = PyDict::new(py);
            d.set_item("homozygous", f.homozygous)?;
            d.set_item("absent_gene", f.absent_gene)?;
            d.set_item("drb_completed", f.drb_completed)?;
            d.set_item("duplicate_read_as_one", f.duplicate_read_as_one)?;
            d.set_item("drb1_conflicts", f.drb1_conflicts)?;
            out.set_item(&f.locus, d)?;
        }
        Ok(out)
    }

    impl Dataset {
        fn per_locus<'py>(
            &self,
            py: Python<'py>,
            counts: Vec<usize>,
        ) -> PyResult<Bound<'py, PyDict>> {
            let d = PyDict::new(py);
            for (locus, n) in self.inner.loci.iter().zip(counts) {
                d.set_item(locus, n)?;
            }
            Ok(d)
        }
    }

    #[pymethods]
    impl Dataset {
        #[getter]
        fn loci(&self) -> Vec<String> {
            self.inner.loci.clone()
        }

        #[getter]
        fn ids(&self) -> Vec<String> {
            self.inner.samples.iter().map(|s| s.id.clone()).collect()
        }

        /// All typed allele names; position in the list is the allele's integer ID.
        #[getter]
        fn alleles(&self) -> Vec<String> {
            self.inner.alleles.names().to_vec()
        }

        /// Lines that could not be read, as (line number, sample id, reason).
        #[getter]
        fn issues(&self) -> Vec<(usize, Option<String>, String)> {
            issues_to_py(&self.inner.issues)
        }

        /// One sample as {locus: GL string or None}.
        fn sample<'py>(&self, py: Python<'py>, index: isize) -> PyResult<Bound<'py, PyDict>> {
            let n = self.inner.samples.len() as isize;
            let i = if index < 0 { index + n } else { index };
            if !(0..n).contains(&i) {
                return Err(PyIndexError::new_err("sample index out of range"));
            }
            let sample = &self.inner.samples[i as usize];
            let d = PyDict::new(py);
            for (locus, obs) in self.inner.loci.iter().zip(&sample.loci) {
                let value = (!obs.is_missing()).then(|| obs.to_gl(&self.inner.alleles));
                d.set_item(locus, value)?;
            }
            Ok(d)
        }

        /// Number of samples with no typing at each locus.
        fn missing_counts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
            self.per_locus(py, self.inner.missing_counts())
        }

        /// Number of samples with only one allele typed at each locus (after filling).
        fn one_allele_counts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
            self.per_locus(py, self.inner.one_allele_counts())
        }

        /// What reading filled in, per locus: {"A": {"homozygous": 12, ...}, ...}.
        fn fill_counts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
            fills_to_py(py, &self.fills)
        }

        fn __len__(&self) -> usize {
            self.inner.samples.len()
        }

        fn __repr__(&self) -> String {
            format!(
                "Dataset({} samples, loci={:?}, {} alleles, {} unreadable lines)",
                self.inner.samples.len(),
                self.inner.loci,
                self.inner.alleles.len(),
                self.inner.issues.len()
            )
        }
    }

    /// Run `read` with the MAC table given as a loaded MacCodes or a path to
    /// alpha.v3.txt; with a path, only the codes used in `path` are loaded.
    fn with_mac<T>(
        path: &PathBuf,
        mac_codes: Option<&Bound<'_, PyAny>>,
        read: impl FnOnce(Option<&hh::MacDictionary>) -> Result<T, hh::Error>,
    ) -> PyResult<T> {
        let loaded;
        let mac: Option<&hh::MacDictionary> = match mac_codes {
            None => None,
            Some(obj) => match obj.cast::<MacCodes>() {
                Ok(m) => Some(&m.get().inner),
                Err(_) => {
                    let mac_path: PathBuf = obj.extract()?;
                    let wanted = hh::collect_mac_codes_in_file(path).map_err(to_py_err)?;
                    if wanted.is_empty() {
                        None // no MAC codes in the input: skip scanning the big table
                    } else {
                        loaded = hh::MacDictionary::from_path_only(mac_path, &wanted)
                            .map_err(to_py_err)?;
                        Some(&loaded)
                    }
                }
            },
        };
        read(mac).map_err(to_py_err)
    }

    /// Shared body of read_mac / read_glsc: read, then fill blanks.
    fn read_with(
        path: PathBuf,
        mac_codes: Option<&Bound<'_, PyAny>>,
        loci: Option<Vec<String>>,
        fill: hh::FillOptions,
        reader: fn(&PathBuf, &hh::ParseOptions) -> Result<hh::Dataset, hh::Error>,
    ) -> PyResult<Dataset> {
        let mut ds = with_mac(&path, mac_codes, |mac| {
            reader(&path, &hh::ParseOptions { loci, mac })
        })?;
        let fills = hh::fill_dataset(&mut ds, &fill);
        Ok(Dataset { inner: ds, fills })
    }

    /// Read a MAC table (tab-separated, two columns per locus).
    ///
    /// `half_typed`: "homozygous" (default) copies a single typed allele across;
    /// "unknown" keeps the second allele unknown. `drb345`: fill DRB3/4/5 from DRB1.
    #[pyfunction]
    #[pyo3(signature = (path, mac_codes=None, loci=None, *, half_typed="homozygous", drb345=true))]
    fn read_mac(
        path: PathBuf,
        mac_codes: Option<&Bound<'_, PyAny>>,
        loci: Option<Vec<String>>,
        half_typed: &str,
        drb345: bool,
    ) -> PyResult<Dataset> {
        let fill = fill_options(half_typed, drb345)?;
        read_with(path, mac_codes, loci, fill, |p, o| hh::read_mac_file(p, o))
    }

    /// Read a GL string table (ID followed by one GL string per locus).
    /// `half_typed` and `drb345` as for read_mac.
    #[pyfunction]
    #[pyo3(signature = (path, mac_codes=None, loci=None, *, half_typed="homozygous", drb345=true))]
    fn read_glsc(
        path: PathBuf,
        mac_codes: Option<&Bound<'_, PyAny>>,
        loci: Option<Vec<String>>,
        half_typed: &str,
        drb345: bool,
    ) -> PyResult<Dataset> {
        let fill = fill_options(half_typed, drb345)?;
        read_with(path, mac_codes, loci, fill, |p, o| hh::read_glsc_file(p, o))
    }

    /// Classify one typing (without locus prefix), e.g. "24:DXTW" -> "mac".
    #[pyfunction]
    fn classify_typing(typing: &str) -> Option<&'static str> {
        hh::classify(typing).map(|k| match k {
            hh::TypingKind::Allele { .. } => "allele",
            hh::TypingKind::GGroup => "G",
            hh::TypingKind::SmallGGroup => "g",
            hh::TypingKind::PGroup => "P",
            hh::TypingKind::Mac => "mac",
        })
    }

    // ------------------------------------------------------------ Nomenclature

    /// IMGT/HLA lookup tables (a data/derived/<release> folder, or Hapl-o-Mat's data/).
    #[pyclass(name = "Nomenclature", module = "hlahaplotools._core", frozen)]
    struct Nomenclature {
        inner: hh::Nomenclature,
    }

    #[pymethods]
    impl Nomenclature {
        #[new]
        fn new(path: PathBuf) -> PyResult<Self> {
            Ok(Self {
                inner: hh::Nomenclature::from_dir(path).map_err(to_py_err)?,
            })
        }

        /// Folder name the tables came from (usually the IMGT/HLA release).
        #[getter]
        fn source(&self) -> Option<String> {
            self.inner.source.clone()
        }

        /// Translate one typing, e.g. ("A*02:01", "g") -> [("A*02:01:01G"..., weight)].
        fn translate(&self, code: &str, resolution: &str) -> PyResult<Vec<(String, f64)>> {
            let res: hh::Resolution = resolution.parse().map_err(PyValueError::new_err)?;
            hh::translate(&self.inner, code, res).map_err(PyValueError::new_err)
        }

        fn __repr__(&self) -> String {
            format!(
                "Nomenclature({:?}: {} g, {} G, {} P groups)",
                self.inner.source.as_deref().unwrap_or("?"),
                self.inner.small_g.len(),
                self.inner.large_g.len(),
                self.inner.p.len()
            )
        }
    }

    // ------------------------------------------------------------ Estimate

    /// Result of estimate().
    #[pyclass(name = "Estimate", module = "hlahaplotools._core", frozen)]
    struct Estimate {
        inner: hh::Estimate,
    }

    #[pymethods]
    impl Estimate {
        #[getter]
        fn loci(&self) -> Vec<String> {
            self.inner.loci.clone()
        }

        #[getter]
        fn resolutions(&self) -> Vec<String> {
            self.inner
                .resolutions
                .iter()
                .map(ToString::to_string)
                .collect()
        }

        /// All haplotypes with nonzero frequency, most frequent first.
        #[getter]
        fn haplotypes(&self) -> Vec<(String, f64)> {
            self.inner.haplotypes.clone()
        }

        #[getter]
        fn iterations(&self) -> usize {
            self.inner.iterations
        }

        #[getter]
        fn converged(&self) -> bool {
            self.inner.converged
        }

        #[getter]
        fn log_likelihood(&self) -> f64 {
            self.inner.log_likelihood
        }

        /// Per iteration: (largest frequency change, log-likelihood).
        #[getter]
        fn history(&self) -> Vec<(f64, f64)> {
            self.inner.history.clone()
        }

        #[getter]
        fn n_samples(&self) -> usize {
            self.inner.n_samples
        }

        /// Samples used that had an untyped or half-typed locus.
        #[getter]
        fn n_samples_incomplete(&self) -> usize {
            self.inner.n_samples_incomplete
        }

        #[getter]
        fn n_phenotypes(&self) -> usize {
            self.inner.n_phenotypes
        }

        #[getter]
        fn n_candidates(&self) -> usize {
            self.inner.n_candidates
        }

        #[getter]
        fn n_completions(&self) -> usize {
            self.inner.n_completions
        }

        /// Seconds spent per stage: {"resolve": s, "model": s, "em": s}.
        #[getter]
        fn timings(&self) -> HashMap<String, f64> {
            self.inner
                .timings
                .iter()
                .map(|&(k, v)| (k.to_string(), v))
                .collect()
        }

        /// Samples left out of the estimate, as (0, sample id, reason).
        #[getter]
        fn issues(&self) -> Vec<(usize, Option<String>, String)> {
            issues_to_py(&self.inner.issues)
        }

        /// Most likely haplotype pair per sample: (sample ID, haplotype 1, haplotype 2, probability).
        fn best_diplotypes(&self) -> Vec<(String, String, String, f64)> {
            self.inner
                .best_diplotypes()
                .iter()
                .map(|d| {
                    let [h1, h2] = d.haplotypes.clone();
                    (d.sample_id.clone(), h1, h2, d.probability)
                })
                .collect()
        }

        /// Haplotypes above `cut`, rescaled to sum to 1 unless renormalize=False.
        #[pyo3(signature = (cut=1e-6, renormalize=true))]
        fn frequencies(&self, cut: f64, renormalize: bool) -> Vec<(String, f64)> {
            self.inner.frequencies(cut, renormalize)
        }

        /// Write `haplotype<TAB>frequency` lines (Hapl-o-Mat's hfs.dat format).
        #[pyo3(signature = (path, cut=1e-6, renormalize=true))]
        fn write_frequencies(&self, path: PathBuf, cut: f64, renormalize: bool) -> PyResult<()> {
            hh::output::write_frequencies(path, &self.inner, cut, renormalize).map_err(to_py_err)
        }

        /// Write the convergence trace (Hapl-o-Mat's epsilon.dat format).
        fn write_history(&self, path: PathBuf) -> PyResult<()> {
            hh::output::write_history(path, &self.inner).map_err(to_py_err)
        }

        fn __repr__(&self) -> String {
            format!(
                "Estimate({} haplotypes from {} samples ({} incomplete), {} iterations, converged={}, logL={:.6})",
                self.inner.haplotypes.len(),
                self.inner.n_samples,
                self.inner.n_samples_incomplete,
                self.inner.iterations,
                self.inner.converged,
                self.inner.log_likelihood
            )
        }
    }

    /// Accept resolutions as a list (one per locus) or a {locus: resolution} dict.
    fn parse_resolutions(loci: &[String], obj: &Bound<'_, PyAny>) -> PyResult<Vec<hh::Resolution>> {
        let names: Vec<String> = if let Ok(map) = obj.extract::<HashMap<String, String>>() {
            loci.iter()
                .map(|l| {
                    map.get(l).cloned().ok_or_else(|| {
                        PyValueError::new_err(format!("no resolution for locus {l}"))
                    })
                })
                .collect::<PyResult<_>>()?
        } else if let Ok(single) = obj.extract::<String>() {
            vec![single; loci.len()]
        } else {
            obj.extract::<Vec<String>>()?
        };
        names
            .iter()
            .map(|s| s.parse().map_err(PyValueError::new_err))
            .collect()
    }

    fn resolutions_or_g(
        loci: &[String],
        obj: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Vec<hh::Resolution>> {
        match obj {
            Some(obj) => parse_resolutions(loci, obj),
            None => Ok(vec![hh::Resolution::SmallG; loci.len()]),
        }
    }

    /// Run `work` without holding the GIL, on `threads` threads if given.
    fn run_detached<T: Send>(
        py: Python<'_>,
        threads: Option<usize>,
        work: impl FnOnce() -> Result<T, hh::Error> + Send,
    ) -> PyResult<T> {
        py.detach(|| match threads {
            Some(n) => rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| hh::Error::Format(e.to_string()))
                .and_then(|pool| pool.install(work)),
            None => work(),
        })
        .map_err(to_py_err)
    }

    #[allow(clippy::too_many_arguments)]
    fn em_options(
        epsilon: f64,
        max_iterations: usize,
        init: &str,
        seed: u64,
        min_genotype_frequency: f64,
        max_completions: usize,
    ) -> PyResult<hh::EmOptions> {
        Ok(hh::EmOptions {
            epsilon,
            max_iterations,
            init: init.parse().map_err(PyValueError::new_err)?,
            seed,
            min_genotype_frequency,
            max_completions,
        })
    }

    /// Estimate haplotype frequencies by EM.
    ///
    /// `resolutions`: "g" for all loci (the default), a list per locus, or {"A": "g", ...}.
    /// `threads`: worker threads (default: all cores). Results do not depend on it.
    #[pyfunction]
    #[pyo3(signature = (
        dataset, nomenclature, resolutions=None, *, epsilon=1e-6, max_iterations=100_000,
        init="numberOccurrence", seed=0, min_genotype_frequency=1e-5, max_completions=10_000, threads=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn estimate(
        py: Python<'_>,
        dataset: &Dataset,
        nomenclature: &Nomenclature,
        resolutions: Option<&Bound<'_, PyAny>>,
        epsilon: f64,
        max_iterations: usize,
        init: &str,
        seed: u64,
        min_genotype_frequency: f64,
        max_completions: usize,
        threads: Option<usize>,
    ) -> PyResult<Estimate> {
        let res = resolutions_or_g(&dataset.inner.loci, resolutions)?;
        let opts = em_options(
            epsilon,
            max_iterations,
            init,
            seed,
            min_genotype_frequency,
            max_completions,
        )?;
        let (ds, nom) = (&dataset.inner, &nomenclature.inner);
        let inner = run_detached(py, threads, || hh::estimate(ds, nom, &res, &opts))?;
        Ok(Estimate { inner })
    }

    // ------------------------------------------------------------ families

    /// A family table: the typing data plus who belongs to which family.
    #[pyclass(name = "FamilyData", module = "hlahaplotools._core", frozen)]
    struct FamilyData {
        ds: hh::Dataset,
        ped: hh::Pedigree,
        fills: Vec<hh::FillCount>,
    }

    #[pymethods]
    impl FamilyData {
        /// The typing data (all members, after filling), as a Dataset.
        #[getter]
        fn dataset(&self) -> Dataset {
            Dataset {
                inner: self.ds.clone(),
                fills: self.fills.clone(),
            }
        }

        /// What reading filled in, per locus (see Dataset.fill_counts).
        fn fill_counts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
            fills_to_py(py, &self.fills)
        }

        #[getter]
        fn loci(&self) -> Vec<String> {
            self.ds.loci.clone()
        }

        /// Families as dicts: {"family", "father", "mother", "children"} (sample IDs).
        fn families<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
            let id = |i: usize| self.ds.samples[i].id.clone();
            self.ped
                .families
                .iter()
                .map(|f| {
                    let d = PyDict::new(py);
                    d.set_item("family", &f.id)?;
                    d.set_item("father", f.father.map(id))?;
                    d.set_item("mother", f.mother.map(id))?;
                    d.set_item(
                        "children",
                        f.children.iter().map(|&c| id(c)).collect::<Vec<_>>(),
                    )?;
                    Ok(d)
                })
                .collect()
        }

        /// Sample IDs not in any family (used as unrelated people).
        #[getter]
        fn unrelated(&self) -> Vec<String> {
            self.ped
                .unrelated
                .iter()
                .map(|&i| self.ds.samples[i].id.clone())
                .collect()
        }

        /// Unreadable rows and members left out of families.
        #[getter]
        fn issues(&self) -> Vec<(usize, Option<String>, String)> {
            let mut all = issues_to_py(&self.ds.issues);
            all.extend(issues_to_py(&self.ped.issues));
            all
        }

        fn __len__(&self) -> usize {
            self.ds.samples.len()
        }

        fn __repr__(&self) -> String {
            format!(
                "FamilyData({} people, {} families, {} unrelated, loci={:?})",
                self.ds.samples.len(),
                self.ped.families.len(),
                self.ped.unrelated.len(),
                self.ds.loci
            )
        }
    }

    /// Read a family table (CSV/TSV with FAMILY_ID, Family_Member and LOCUS_1/LOCUS_2 columns).
    ///
    /// `mac_codes`: a MacCodes or a path to alpha.v3.txt, if the table uses MAC codes.
    /// `half_typed` and `drb345` as for read_mac.
    #[pyfunction]
    #[pyo3(signature = (path, mac_codes=None, loci=None, *, half_typed="homozygous", drb345=true))]
    fn read_family_table(
        path: PathBuf,
        mac_codes: Option<&Bound<'_, PyAny>>,
        loci: Option<Vec<String>>,
        half_typed: &str,
        drb345: bool,
    ) -> PyResult<FamilyData> {
        let fill = fill_options(half_typed, drb345)?;
        let (mut ds, ped) = with_mac(&path, mac_codes, |mac| {
            hh::read_family_table(&path, &hh::ParseOptions { loci, mac })
        })?;
        let fills = hh::fill_dataset(&mut ds, &fill);
        Ok(FamilyData { ds, ped, fills })
    }

    /// Check each child against its parents at the target resolution.
    /// Returns dicts {"family", "sample", "locus", "message"}.
    #[pyfunction]
    #[pyo3(signature = (data, nomenclature, resolutions=None))]
    fn check_mendelian<'py>(
        py: Python<'py>,
        data: &FamilyData,
        nomenclature: &Nomenclature,
        resolutions: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let res = resolutions_or_g(&data.ds.loci, resolutions)?;
        let issues = hh::check_mendelian(&data.ds, &data.ped, &nomenclature.inner, &res)
            .map_err(to_py_err)?;
        issues
            .into_iter()
            .map(|i| {
                let d = PyDict::new(py);
                d.set_item("family", i.family)?;
                d.set_item("sample", i.sample_id)?;
                d.set_item("locus", i.locus)?;
                d.set_item("message", i.message)?;
                Ok(d)
            })
            .collect()
    }

    /// Check every family for members whose typing does not fit their stated
    /// relationships. Returns dicts {"family", "sample", "relative", "kind",
    /// "loci", "message"}; kinds: structure, parent_excluded,
    /// single_locus_mismatch, trio_inconsistent, possible_duplicate,
    /// not_full_sibling, recombination, inheritance_mismatch.
    #[pyfunction]
    #[pyo3(signature = (
        data, nomenclature, resolutions=None, *, min_exclusions=2, duplicate_fraction=0.9, min_loci_compared=5,
        cross_family_duplicates=false
    ))]
    #[allow(clippy::too_many_arguments)]
    fn audit_families<'py>(
        py: Python<'py>,
        data: &FamilyData,
        nomenclature: &Nomenclature,
        resolutions: Option<&Bound<'_, PyAny>>,
        min_exclusions: usize,
        duplicate_fraction: f64,
        min_loci_compared: usize,
        cross_family_duplicates: bool,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let res = resolutions_or_g(&data.ds.loci, resolutions)?;
        let opts = hh::AuditOptions {
            min_exclusions,
            duplicate_fraction,
            min_loci_compared,
            cross_family_duplicates,
        };
        let flags = hh::audit_families(&data.ds, &data.ped, &nomenclature.inner, &res, &opts)
            .map_err(to_py_err)?;
        flags
            .into_iter()
            .map(|f| {
                let d = PyDict::new(py);
                d.set_item("family", f.family)?;
                d.set_item("sample", f.sample_id)?;
                d.set_item("relative", f.relative)?;
                d.set_item("kind", f.kind.as_str())?;
                d.set_item("loci", f.loci)?;
                d.set_item("message", f.message)?;
                Ok(d)
            })
            .collect()
    }

    /// Phasing of one family: the parents' four haplotypes and each child's inheritance.
    #[pyclass(name = "FamilyPhase", module = "hlahaplotools._core", frozen)]
    struct FamilyPhase {
        inner: hh::FamilyPhase,
    }

    #[pymethods]
    impl FamilyPhase {
        #[getter]
        fn family(&self) -> String {
            self.inner.family.clone()
        }
        #[getter]
        fn father(&self) -> Option<String> {
            self.inner.father.clone()
        }
        #[getter]
        fn mother(&self) -> Option<String> {
            self.inner.mother.clone()
        }
        /// Haplotypes [A, B, C, D]: A and B are the father's, C and D the mother's;
        /// A and C are those of the lowest-numbered typed child. Unknown loci appear as "LOCUS*?".
        #[getter]
        fn haplotypes(&self) -> Vec<String> {
            self.inner.haplotypes.to_vec()
        }
        /// The same haplotypes keyed by label: {"A": ..., "B": ..., "C": ..., "D": ...}.
        #[getter]
        fn labelled_haplotypes(&self) -> HashMap<String, String> {
            ["A", "B", "C", "D"]
                .iter()
                .zip(&self.inner.haplotypes)
                .map(|(l, h)| (l.to_string(), h.clone()))
                .collect()
        }
        /// (child sample ID, "A"/"B" from the father, "C"/"D" from the mother),
        /// in child-number order; None where the data do not decide.
        #[getter]
        fn children(&self) -> Vec<(String, Option<String>, Option<String>)> {
            self.inner
                .children
                .iter()
                .map(|c| {
                    (
                        c.sample_id.clone(),
                        c.paternal.map(String::from),
                        c.maternal.map(String::from),
                    )
                })
                .collect()
        }
        #[getter]
        fn inconsistent_loci(&self) -> Vec<String> {
            self.inner.inconsistent_loci.clone()
        }
        #[getter]
        fn ambiguous_loci(&self) -> Vec<String> {
            self.inner.ambiguous_loci.clone()
        }
        /// "complete", "partial" or "failed".
        #[getter]
        fn status(&self) -> &'static str {
            match self.inner.status {
                hh::PhaseStatus::Complete => "complete",
                hh::PhaseStatus::Partial => "partial",
                hh::PhaseStatus::Failed(_) => "failed",
            }
        }
        #[getter]
        fn notes(&self) -> Vec<String> {
            let mut notes = self.inner.notes.clone();
            if let hh::PhaseStatus::Failed(reason) = &self.inner.status {
                notes.insert(0, reason.clone());
            }
            notes
        }
        /// ID used for a parent in estimates: the typed ID or "family/father (inferred)".
        fn parent_id(&self, parent: usize) -> PyResult<String> {
            if parent > 1 {
                return Err(PyValueError::new_err(
                    "parent must be 0 (father) or 1 (mother)",
                ));
            }
            Ok(self.inner.parent_id(parent))
        }
        fn __repr__(&self) -> String {
            format!("FamilyPhase({:?}, {})", self.inner.family, self.status())
        }
    }

    /// Phase every family at the target resolution.
    #[pyfunction]
    #[pyo3(signature = (data, nomenclature, resolutions=None))]
    fn phase_families(
        data: &FamilyData,
        nomenclature: &Nomenclature,
        resolutions: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Vec<FamilyPhase>> {
        let res = resolutions_or_g(&data.ds.loci, resolutions)?;
        let phases = hh::phase_families(&data.ds, &data.ped, &nomenclature.inner, &res)
            .map_err(to_py_err)?;
        Ok(phases
            .into_iter()
            .map(|inner| FamilyPhase { inner })
            .collect())
    }

    /// Estimate from founders (parents) only; children are used to phase them.
    ///
    /// With `use_family_phase=False`, founders enter as unrelated people
    /// (useful for compare_to_families). Returns (Estimate, [FamilyPhase]).
    #[pyfunction]
    #[pyo3(signature = (
        data, nomenclature, resolutions=None, *, use_family_phase=true, epsilon=1e-6, max_iterations=100_000,
        init="numberOccurrence", seed=0, min_genotype_frequency=1e-5, max_completions=10_000, threads=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn estimate_families(
        py: Python<'_>,
        data: &FamilyData,
        nomenclature: &Nomenclature,
        resolutions: Option<&Bound<'_, PyAny>>,
        use_family_phase: bool,
        epsilon: f64,
        max_iterations: usize,
        init: &str,
        seed: u64,
        min_genotype_frequency: f64,
        max_completions: usize,
        threads: Option<usize>,
    ) -> PyResult<(Estimate, Vec<FamilyPhase>)> {
        let res = resolutions_or_g(&data.ds.loci, resolutions)?;
        let opts = em_options(
            epsilon,
            max_iterations,
            init,
            seed,
            min_genotype_frequency,
            max_completions,
        )?;
        let (ds, ped, nom) = (&data.ds, &data.ped, &nomenclature.inner);
        let (est, phases) = run_detached(py, threads, || {
            hh::estimate_families(ds, ped, nom, &res, &opts, use_family_phase)
        })?;
        Ok((
            Estimate { inner: est },
            phases
                .into_iter()
                .map(|inner| FamilyPhase { inner })
                .collect(),
        ))
    }

    /// Compare an estimate with family phasing. Returns a dict with "founders"
    /// (per-founder comparison), "n_exact", "mean_locus_agreement",
    /// "frequencies" [(haplotype, family freq, EM freq)], "total_variation"
    /// and "n_family_haplotypes".
    #[pyfunction]
    fn compare_to_families<'py>(
        py: Python<'py>,
        estimate: &Estimate,
        phases: Vec<PyRef<'_, FamilyPhase>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let phases: Vec<hh::FamilyPhase> = phases.iter().map(|p| p.inner.clone()).collect();
        let cmp = hh::compare_to_families(&estimate.inner, &phases);
        let founders: Vec<Bound<'py, PyDict>> = cmp
            .founders
            .iter()
            .map(|f| {
                let d = PyDict::new(py);
                d.set_item("family", &f.family)?;
                d.set_item("sample", &f.sample_id)?;
                d.set_item("family_haplotypes", f.family_haplotypes.to_vec())?;
                d.set_item("em_haplotypes", f.em_haplotypes.to_vec())?;
                d.set_item("em_probability", f.em_probability)?;
                d.set_item("exact", f.exact)?;
                d.set_item("loci_agreeing", f.loci_agreeing)?;
                d.set_item("n_loci", f.n_loci)?;
                Ok(d)
            })
            .collect::<PyResult<_>>()?;
        let out = PyDict::new(py);
        out.set_item("founders", founders)?;
        out.set_item("n_compared", cmp.founders.len())?;
        out.set_item("n_exact", cmp.n_exact)?;
        out.set_item("mean_locus_agreement", cmp.mean_locus_agreement)?;
        out.set_item("frequencies", cmp.frequencies)?;
        out.set_item("total_variation", cmp.total_variation)?;
        out.set_item("n_family_haplotypes", cmp.n_family_haplotypes)?;
        Ok(out)
    }
}
