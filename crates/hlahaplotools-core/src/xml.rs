//! Readers for XML typing reports: NMDP HML and GenDx NGSengine exports.
//!
//! Both are read as a stream, so file size does not limit memory. Each locus
//! result becomes a GL string and goes through the same path as GLSC input,
//! keeping all ambiguity (`|`, `/`) for the EM to weigh.
//!
//! **HML** (1.0 and 1.0.1, namespace not required): every `<sample id=...>`,
//! with one GL string per `<allele-assignment>`, from `<glstring>` or, when
//! there is none, from its `<haploid locus=... type=...>` elements.
//!
//! **NGSengine** (`<ProjectXml>`, the "TARR" export): every `<Locus>` with its
//! `<SampleName>`, the reviewed result in
//! `<ConclusiveTypingResult>/<Typing>/<GLString>`, and the review status in
//! `<ReviewList CurrentApprovalStatus=...>`. Loci not approved (e.g.
//! `Rejected`) are left untyped and listed in [`Dataset::notes`], unless
//! [`XmlOptions::include_unapproved`] is set. Loci with no review list are read.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use quick_xml::Reader;
use quick_xml::events::{BytesRef, BytesStart, Event};

use crate::Error;
use crate::dataset::{Dataset, Issue};
use crate::parse::{GlRecord, ParseOptions, dataset_from_gl};

/// Which XML layout a file uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XmlKind {
    Hml,
    NgsEngine,
}

impl std::str::FromStr for XmlKind {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "hml" => Ok(XmlKind::Hml),
            "ngsengine" | "tarr" => Ok(XmlKind::NgsEngine),
            other => Err(format!(
                "unknown XML format {other:?} (use hml or ngsengine)"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct XmlOptions {
    /// NGSengine: also read loci whose review status is not `Approved`.
    pub include_unapproved: bool,
}

/// Guess the layout from the start of a file: `<hml` or `<ProjectXml`.
pub fn sniff_xml(path: impl AsRef<Path>) -> Result<Option<XmlKind>, Error> {
    use std::io::Read;
    let mut head = vec![0u8; 8192];
    let n = File::open(path)?.read(&mut head)?;
    let head = String::from_utf8_lossy(&head[..n]);
    Ok(
        if head.contains("<ProjectXml") || head.contains("NGSengine") {
            Some(XmlKind::NgsEngine)
        } else if head.contains("<hml") || head.contains(":hml") {
            Some(XmlKind::Hml)
        } else {
            None
        },
    )
}

/// Read one or more XML files (e.g. one HML file per sample) into one dataset.
pub fn read_xml_files<P: AsRef<Path>>(
    paths: &[P],
    kind: XmlKind,
    xml: &XmlOptions,
    opts: &ParseOptions,
) -> Result<Dataset, Error> {
    let mut collected = Collected::default();
    for path in paths {
        let path = path.as_ref();
        let reader = BufReader::with_capacity(1 << 20, File::open(path)?);
        let label = path.display().to_string();
        match kind {
            XmlKind::Hml => collect_hml(reader, &label, &mut collected)?,
            XmlKind::NgsEngine => collect_ngsengine(reader, &label, xml, &mut collected)?,
        }
    }
    Ok(collected.into_dataset(opts))
}

/// Read one XML document from a reader.
pub fn read_xml<R: BufRead>(
    reader: R,
    kind: XmlKind,
    xml: &XmlOptions,
    opts: &ParseOptions,
) -> Result<Dataset, Error> {
    let mut collected = Collected::default();
    match kind {
        XmlKind::Hml => collect_hml(reader, "input", &mut collected)?,
        XmlKind::NgsEngine => collect_ngsengine(reader, "input", xml, &mut collected)?,
    }
    Ok(collected.into_dataset(opts))
}

// ---------------------------------------------------------------- shared

/// Samples gathered across files, in first-seen order.
#[derive(Default)]
struct Collected {
    order: Vec<String>,
    /// Sample -> (locus -> GL strings seen for it).
    loci: HashMap<String, Vec<(String, String)>>,
    issues: Vec<Issue>,
    notes: Vec<Issue>,
}

impl Collected {
    /// Add a GL string for a sample; `^`-joined loci are kept apart.
    fn add(&mut self, sample: &str, gl: String) {
        if !self.loci.contains_key(sample) {
            self.order.push(sample.to_owned());
        }
        let entry = self.loci.entry(sample.to_owned()).or_default();
        let gl: String = gl
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>()
            .replace("HLA-", "");
        for part in gl.split('^').filter(|p| !p.is_empty()) {
            let locus = part
                .split_once('*')
                .map_or_else(String::new, |(l, _)| l.to_owned());
            entry.push((locus, part.to_owned()));
        }
    }

    fn note(&mut self, sample: &str, message: String) {
        self.notes.push(Issue {
            line: 0,
            sample_id: Some(sample.to_owned()),
            message,
        });
    }

    fn into_dataset(mut self, opts: &ParseOptions) -> Dataset {
        let mut records = Vec::with_capacity(self.order.len());
        for id in std::mem::take(&mut self.order) {
            let entries = self.loci.remove(&id).unwrap_or_default();
            // The same locus reported twice: keep it once if the results agree,
            // otherwise leave it untyped (we cannot tell which is right).
            let mut by_locus: Vec<(String, Vec<String>)> = Vec::new();
            for (locus, compact) in entries {
                match by_locus
                    .iter_mut()
                    .find(|(l, _)| *l == locus && !l.is_empty())
                {
                    Some((_, gls)) => {
                        if !gls.contains(&compact) {
                            gls.push(compact);
                        }
                    }
                    None => by_locus.push((locus, vec![compact])),
                }
            }
            let mut fields = Vec::new();
            for (locus, gls) in by_locus {
                if gls.len() == 1 {
                    fields.extend(gls);
                } else {
                    self.notes.push(Issue {
                        line: 0,
                        sample_id: Some(id.clone()),
                        message: format!(
                            "{locus}: {} different results ({}); left untyped",
                            gls.len(),
                            gls.join(" vs ")
                        ),
                    });
                }
            }
            records.push(GlRecord {
                line: 0,
                id,
                fields,
            });
        }
        dataset_from_gl(records, self.issues, self.notes, opts)
    }
}

fn xml_error(label: &str, reader_pos: u64, e: impl std::fmt::Display) -> Error {
    Error::Format(format!("{label}: XML error at byte {reader_pos}: {e}"))
}

fn attr(e: &BytesStart, name: &str) -> Option<String> {
    e.attributes()
        .flatten()
        .find(|a| a.key.local_name().as_ref() == name)
        .map(|a| {
            a.normalized_value(quick_xml::XmlVersion::Implicit1_0)
                .map(|v| v.into_owned())
                .unwrap_or_else(|_| a.value.clone().into_owned())
        })
}

/// `&lt;`, `&#65;` and similar, as text.
fn entity(r: &BytesRef) -> String {
    if let Ok(Some(c)) = r.resolve_char_ref() {
        return c.to_string();
    }
    match r.as_ref() {
        "lt" => "<",
        "gt" => ">",
        "amp" => "&",
        "apos" => "'",
        "quot" => "\"",
        _ => "",
    }
    .to_owned()
}

fn local(e: &BytesStart) -> String {
    e.local_name().as_ref().to_owned()
}

// ---------------------------------------------------------------- HML

fn collect_hml<R: BufRead>(reader: R, label: &str, out: &mut Collected) -> Result<(), Error> {
    let mut xml = Reader::from_reader(reader);
    let mut buf = Vec::new();
    let mut stack: Vec<String> = Vec::new();

    let mut sample: Option<String> = None;
    let mut in_glstring = false;
    let mut text = String::new();
    // Per allele-assignment: GL strings and haploid alleles found in it.
    let mut assignment_gls: Vec<String> = Vec::new();
    let mut haploids: Vec<(String, String)> = Vec::new();

    loop {
        let pos = xml.buffer_position();
        let event = xml
            .read_event_into(&mut buf)
            .map_err(|e| xml_error(label, pos, e))?;
        match event {
            Event::Start(ref e) | Event::Empty(ref e) => {
                let name = local(e);
                let empty = matches!(event, Event::Empty(_));
                match name.as_str() {
                    "sample" => {
                        let id = attr(e, "id").unwrap_or_default();
                        if id.is_empty() {
                            out.issues.push(Issue {
                                line: 0,
                                sample_id: None,
                                message: format!("{label}: <sample> without an id attribute"),
                            });
                        }
                        sample = Some(id);
                    }
                    "allele-assignment" => {
                        assignment_gls.clear();
                        haploids.clear();
                    }
                    "glstring" => {
                        text.clear();
                        in_glstring = !empty;
                        if empty
                            && let Some(uri) = attr(e, "uri")
                            && let Some(s) = &sample
                        {
                            out.note(
                                s,
                                format!(
                                    "GL string given only as a URI ({uri}); not read. Resolve it \
                                     with the GL service and give the GL string as text"
                                ),
                            );
                        }
                    }
                    "haploid" => {
                        if let (Some(locus), Some(ty)) = (attr(e, "locus"), attr(e, "type")) {
                            haploids.push((locus, ty));
                        }
                    }
                    _ => {}
                }
                if !empty {
                    stack.push(name);
                }
            }
            Event::Text(t) if in_glstring => text.push_str(t.as_ref()),
            Event::CData(t) if in_glstring => text.push_str(&t),
            Event::GeneralRef(r) if in_glstring => text.push_str(&entity(&r)),
            Event::End(e) => {
                let name = e.local_name().as_ref().to_owned();
                stack.pop();
                match name.as_str() {
                    "glstring" => {
                        in_glstring = false;
                        if !text.trim().is_empty() {
                            assignment_gls.push(std::mem::take(&mut text));
                        }
                    }
                    "allele-assignment" => {
                        if let Some(s) = sample.clone().filter(|s| !s.is_empty()) {
                            if !assignment_gls.is_empty() {
                                for gl in assignment_gls.drain(..) {
                                    out.add(&s, gl);
                                }
                            } else if !haploids.is_empty() {
                                out.add(&s, haploid_gl(&haploids));
                            }
                        }
                        assignment_gls.clear();
                        haploids.clear();
                    }
                    "sample" => {
                        if let Some(s) = sample.take().filter(|s| !s.is_empty())
                            && !out.loci.contains_key(&s)
                        {
                            out.note(&s, "no typing results in the HML file".into());
                        }
                    }
                    _ => {}
                }
            }
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(())
}

/// `<haploid locus="HLA-A" type="01:01"/>` pairs -> `A*01:01+A*02:01`, per locus.
fn haploid_gl(haploids: &[(String, String)]) -> String {
    let mut per_locus: Vec<(String, Vec<String>)> = Vec::new();
    for (locus, ty) in haploids {
        let locus = locus.replace("HLA-", "");
        let allele = if ty.contains('*') {
            ty.replace("HLA-", "")
        } else {
            format!("{locus}*{ty}")
        };
        match per_locus.iter_mut().find(|(l, _)| *l == locus) {
            Some((_, a)) => a.push(allele),
            None => per_locus.push((locus, vec![allele])),
        }
    }
    per_locus
        .into_iter()
        .map(|(_, a)| match a.len() {
            1 => format!("{}+", a[0]),
            _ => format!("{}+{}", a[0], a[1]),
        })
        .collect::<Vec<_>>()
        .join("^")
}

// ---------------------------------------------------------------- NGSengine

#[derive(Default)]
struct LocusResult {
    gene: String,
    sample: String,
    status: Option<String>,
    gl: Option<String>,
}

fn collect_ngsengine<R: BufRead>(
    reader: R,
    label: &str,
    opts: &XmlOptions,
    out: &mut Collected,
) -> Result<(), Error> {
    let mut xml = Reader::from_reader(reader);
    let mut buf = Vec::new();
    let mut stack: Vec<String> = Vec::new();

    let mut sample_name = String::new();
    let mut locus: Option<LocusResult> = None;
    let mut text = String::new();

    // Text is only kept for the few elements we read.
    let wanted = |stack: &[String]| -> bool {
        let n = stack.len();
        let ends = |tail: &[&str]| {
            n >= tail.len()
                && stack[n - tail.len()..]
                    .iter()
                    .zip(tail)
                    .all(|(a, b)| a == b)
        };
        ends(&["Sample", "Name"])
            || ends(&["Locus", "Name"])
            || ends(&["Locus", "SampleName"])
            || ends(&["ConclusiveTypingResult", "Typing", "GLString"])
    };

    loop {
        let pos = xml.buffer_position();
        let event = xml
            .read_event_into(&mut buf)
            .map_err(|e| xml_error(label, pos, e))?;
        match event {
            Event::Start(ref e) | Event::Empty(ref e) => {
                let name = local(e);
                let empty = matches!(event, Event::Empty(_));
                match name.as_str() {
                    "Locus" if !empty => {
                        locus = Some(LocusResult {
                            sample: sample_name.clone(),
                            ..Default::default()
                        });
                    }
                    "ReviewList" if stack.last().is_some_and(|s| s == "Locus") => {
                        if let Some(l) = locus.as_mut() {
                            l.status = attr(e, "CurrentApprovalStatus");
                        }
                    }
                    _ => {}
                }
                if !empty {
                    stack.push(name);
                    text.clear();
                }
            }
            Event::Text(t) if wanted(&stack) => text.push_str(t.as_ref()),
            Event::GeneralRef(r) if wanted(&stack) => text.push_str(&entity(&r)),
            Event::End(_) => {
                let n = stack.len();
                let is = |tail: &[&str]| {
                    n >= tail.len()
                        && stack[n - tail.len()..]
                            .iter()
                            .zip(tail)
                            .all(|(a, b)| a == b)
                };
                let value = text.trim().to_owned();
                if is(&["Sample", "Name"]) {
                    sample_name = value;
                } else if is(&["Locus", "Name"]) {
                    if let Some(l) = locus.as_mut() {
                        l.gene = value.replace("HLA-", "");
                    }
                } else if is(&["Locus", "SampleName"]) {
                    if let Some(l) = locus.as_mut().filter(|_| !value.is_empty()) {
                        l.sample = value;
                    }
                } else if is(&["ConclusiveTypingResult", "Typing", "GLString"]) {
                    if let Some(l) = locus.as_mut().filter(|_| !value.is_empty()) {
                        l.gl = Some(value);
                    }
                } else if is(&["Locus"])
                    && let Some(l) = locus.take()
                {
                    finish_locus(l, opts, out);
                }
                stack.pop();
                text.clear();
            }
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(())
}

fn finish_locus(l: LocusResult, opts: &XmlOptions, out: &mut Collected) {
    if l.sample.is_empty() {
        out.issues.push(Issue {
            line: 0,
            sample_id: None,
            message: format!("{} result without a sample name", l.gene),
        });
        return;
    }
    let approved = match l.status.as_deref() {
        None => true,
        Some(s) => s.eq_ignore_ascii_case("Approved"),
    };
    // Register the sample even if this locus is not used.
    if !out.loci.contains_key(&l.sample) {
        out.order.push(l.sample.clone());
        out.loci.insert(l.sample.clone(), Vec::new());
    }
    match (&l.gl, approved || opts.include_unapproved) {
        (None, _) => out.note(
            &l.sample,
            format!("{}: no conclusive typing result; left untyped", l.gene),
        ),
        (Some(_), false) => out.note(
            &l.sample,
            format!(
                "{}: result not approved ({}); left untyped",
                l.gene,
                l.status.as_deref().unwrap_or("?")
            ),
        ),
        (Some(gl), true) => out.add(&l.sample, gl.clone()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gl(ds: &Dataset, sample: usize, locus: &str) -> String {
        let k = ds.loci.iter().position(|l| l == locus).unwrap();
        ds.samples[sample].loci[k].to_gl(&ds.alleles)
    }

    const HML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<hml xmlns="http://schemas.nmdp.org/spec/hml/1.0.1" version="1.0.1">
  <sample id="S1">
    <typing gene-family="HLA">
      <allele-assignment allele-db="IMGT/HLA" allele-version="3.53.0">
        <glstring>
          HLA-A*02:01:01G+HLA-A*03:01:01G
        </glstring>
      </allele-assignment>
    </typing>
    <typing gene-family="HLA">
      <allele-assignment>
        <glstring>HLA-B*07:02+HLA-B*44:02|HLA-B*07:02+HLA-B*44:03</glstring>
      </allele-assignment>
    </typing>
    <typing gene-family="HLA">
      <allele-assignment>
        <glstring>HLA-DRB1*15:01+HLA-DRB1*03:01^HLA-DRB5*01:01+HLA-DRB5*NNNN</glstring>
      </allele-assignment>
    </typing>
  </sample>
  <sample id="S2">
    <typing gene-family="HLA">
      <allele-assignment>
        <haploid locus="HLA-A" method="DNA" type="01:01"/>
        <haploid locus="HLA-A" method="DNA" type="24:02"/>
      </allele-assignment>
    </typing>
  </sample>
</hml>"#;

    #[test]
    fn hml_keeps_ambiguity_and_reads_haploids() {
        let ds = read_xml(
            HML.as_bytes(),
            XmlKind::Hml,
            &XmlOptions::default(),
            &ParseOptions::default(),
        )
        .unwrap();
        assert_eq!(ds.loci, ["A", "B", "DRB1", "DRB5"]);
        assert_eq!(ds.samples.len(), 2);
        assert_eq!(ds.samples[0].id, "S1");
        assert_eq!(gl(&ds, 0, "A"), "A*02:01:01G+A*03:01:01G");
        assert_eq!(gl(&ds, 0, "B"), "B*07:02+B*44:02|B*07:02+B*44:03");
        assert_eq!(gl(&ds, 0, "DRB1"), "DRB1*15:01+DRB1*03:01");
        assert_eq!(gl(&ds, 0, "DRB5"), "DRB5*01:01+DRB5*NNNN");
        assert_eq!(gl(&ds, 1, "A"), "A*01:01+A*24:02");
        assert!(ds.samples[1].loci[1].is_missing());
        assert!(ds.issues.is_empty());
    }

    const TARR: &str = r#"<?xml version="1.0" encoding="utf-8"?>
<ProjectXml>
  <Samples>
    <Sample>
      <Name>P1_barcode01</Name>
      <Loci>
        <Locus>
          <Name>HLA-A</Name>
          <SampleName>P1_barcode01</SampleName>
          <TypingResultOptimized><GLString>A*99:99+A*99:99</GLString></TypingResultOptimized>
          <ConclusiveTypingResult>
            <Typing><GLString>A*01:01:01:01+A*02:01:01:01</GLString></Typing>
          </ConclusiveTypingResult>
          <ReviewList CurrentApprovalStatus="Approved" Level="Second reviewer"/>
        </Locus>
        <Locus>
          <Name>HLA-B</Name>
          <SampleName>P1_barcode01</SampleName>
          <ConclusiveTypingResult>
            <Typing><GLString>B*08:01:01:01+B*40:01:02</GLString></Typing>
          </ConclusiveTypingResult>
          <ReviewList CurrentApprovalStatus="Rejected"/>
        </Locus>
        <Locus>
          <Name>HLA-C</Name>
          <SampleName>P1_barcode01</SampleName>
          <ConclusiveTypingResult>
            <Typing><GLString>C*07:01:01:01+C*03:04:01:01/C*03:04:01:02</GLString></Typing>
          </ConclusiveTypingResult>
        </Locus>
      </Loci>
    </Sample>
  </Samples>
</ProjectXml>"#;

    #[test]
    fn ngsengine_reads_conclusive_results_and_skips_rejected() {
        let ds = read_xml(
            TARR.as_bytes(),
            XmlKind::NgsEngine,
            &XmlOptions::default(),
            &ParseOptions::default(),
        )
        .unwrap();
        assert_eq!(ds.samples.len(), 1);
        assert_eq!(ds.samples[0].id, "P1_barcode01");
        assert_eq!(ds.loci, ["A", "C"]);
        assert_eq!(gl(&ds, 0, "A"), "A*01:01:01:01+A*02:01:01:01");
        assert_eq!(gl(&ds, 0, "C"), "C*07:01:01:01+C*03:04:01:01/C*03:04:01:02");
        assert_eq!(ds.notes.len(), 1);
        assert!(ds.notes[0].message.contains("not approved (Rejected)"));

        let all = read_xml(
            TARR.as_bytes(),
            XmlKind::NgsEngine,
            &XmlOptions {
                include_unapproved: true,
            },
            &ParseOptions::default(),
        )
        .unwrap();
        assert_eq!(all.loci, ["A", "B", "C"]);
        assert_eq!(gl(&all, 0, "B"), "B*08:01:01:01+B*40:01:02");
    }

    #[test]
    fn conflicting_repeats_are_left_untyped() {
        let doc = TARR
            .replace("B*08:01:01:01+B*40:01:02", "A*01:01:01:01+A*03:01:01:01")
            .replace("HLA-B", "HLA-A");
        let ds = read_xml(
            doc.as_bytes(),
            XmlKind::NgsEngine,
            &XmlOptions {
                include_unapproved: true,
            },
            &ParseOptions {
                loci: Some(vec!["A".into(), "C".into()]),
                mac: None,
            },
        )
        .unwrap();
        assert!(ds.samples[0].loci[0].is_missing());
        assert!(!ds.samples[0].loci[1].is_missing());
        assert!(
            ds.notes
                .iter()
                .any(|n| n.message.contains("2 different results"))
        );
    }
}
