"""Reading GLS (pull + GLID), HML and NGSengine XML input. Fixtures are synthetic."""

import json
import subprocess
import sys

from conftest import DERIVED, needs_derived

import hlahaplotools as h

GLID = "1;A*01:01+A*02:01\n2;B*08:01+B*07:02\n3;HLA-A*01:01+HLA-A*01:01\n4;B*08:01+B*08:01\n"
PULL = "p1;1:2\np2;3:4\np3;3:0\n"

HML = """<?xml version="1.0" encoding="UTF-8"?>
<hml xmlns="http://schemas.nmdp.org/spec/hml/1.0.1" version="1.0.1">
  <sample id="{id}">
    <typing gene-family="HLA"><allele-assignment allele-db="IMGT/HLA" allele-version="3.53.0">
      <glstring>
        HLA-A*{a}
      </glstring>
    </allele-assignment></typing>
    <typing gene-family="HLA"><allele-assignment><glstring>HLA-B*{b}</glstring></allele-assignment></typing>
  </sample>
</hml>
"""

TARR = """<?xml version="1.0" encoding="utf-8"?>
<ProjectXml xmlns:d1p1="http://www.w3.org/2001/XMLSchema-instance">
  <AnalysisSoftware><Software>NGSengine</Software></AnalysisSoftware>
  <Samples>{samples}</Samples>
</ProjectXml>
"""
LOCUS = """<Locus><Name>HLA-{gene}</Name><SampleName>{id}</SampleName>
  <ConclusiveTypingResult><Typing><GLString>{gl}</GLString></Typing></ConclusiveTypingResult>
  <ReviewList CurrentApprovalStatus="{status}" Level="Second reviewer"/></Locus>"""


def tarr(rows):
    samples = []
    for sid, loci in rows:
        body = "".join(LOCUS.format(gene=g, id=sid, gl=gl, status=st) for g, gl, st in loci)
        samples.append(f"<Sample><Name>{sid}</Name><Loci>{body}</Loci></Sample>")
    return TARR.format(samples="".join(samples))


def test_read_gls_keeps_untyped_loci(tmp_path):
    (tmp_path / "x.glid").write_text(GLID)
    (tmp_path / "x.pull").write_text(PULL)
    ds = h.read_gls(tmp_path / "x.pull", tmp_path / "x.glid")
    assert ds.loci == ["A", "B"] and len(ds) == 3
    assert ds.sample(1) == {"A": "A*01:01+A*01:01", "B": "B*08:01+B*08:01"}
    assert ds.sample(2)["B"] is None  # GL ID 0: untyped, not dropped
    assert ds.missing_counts() == {"A": 0, "B": 1}


def test_read_hml_files_and_folder(tmp_path):
    d = tmp_path / "hml"
    d.mkdir()
    (d / "s1.hml").write_text(HML.format(id="s1", a="01:01+HLA-A*02:01|HLA-A*01:01+HLA-A*02:05", b="08:01+HLA-B*07:02"))
    (d / "s2.hml").write_text(HML.format(id="s2", a="03:01+HLA-A*11:01", b="35:01"))
    ds = h.read_hml(sorted(d.iterdir()))
    assert ds.ids == ["s1", "s2"]
    assert ds.sample(0)["A"] == "A*01:01+A*02:01|A*01:01+A*02:05"
    # One B allele reported, no '+': read as one allele typed, then homozygous.
    assert ds.sample(1)["B"] == "B*35:01+B*35:01"
    assert h.sniff_xml(d / "s1.hml") == "hml"

    run = subprocess.run([sys.executable, "-m", "hlahaplotools", "summary", str(d)], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert "Dataset(2 samples" in run.stdout and "(HML)" in run.stdout


def test_read_ngsengine_uses_reviewed_results(tmp_path):
    p = tmp_path / "run.xml"
    p.write_text(
        tarr(
            [
                ("bc01", [("A", "A*01:01:01:01+A*02:01:01:01", "Approved"), ("B", "B*08:01+B*07:02", "Rejected")]),
                ("bc02", [("A", "A*03:01:01:01+A*03:01:01:01", "Approved"), ("B", "B*44:02+B*07:02", "Approved")]),
            ]
        )
    )
    assert h.sniff_xml(p) == "ngsengine"
    ds = h.read_ngsengine(p)
    assert ds.ids == ["bc01", "bc02"]
    assert ds.sample(0)["B"] is None
    assert [(sid, "Rejected" in why) for _, sid, why in ds.notes] == [("bc01", True)]
    every = h.read_ngsengine(p, include_unapproved=True)
    assert every.sample(0)["B"] == "B*08:01+B*07:02" and not every.notes


@needs_derived
def test_cli_estimate_from_ngsengine(tmp_path):
    rows = []
    for i in range(30):
        a, b = ("01:01", "08:01") if i % 3 else ("02:01", "07:02")
        rows.append((f"bc{i:02}", [("A", f"A*{a}+A*{a}", "Approved"), ("B", f"B*{b}+B*{b}", "Approved")]))
    rows.append(("bc99", [("A", "A*01:01+A*02:01", "Rejected"), ("B", "B*08:01+B*07:02", "Approved")]))
    p = tmp_path / "run.tarr.xml"
    p.write_text(tarr(rows))
    out = tmp_path / "res"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "hlahaplotools",
            "estimate",
            str(p),
            "--loci",
            "A:2f,B:2f",
            "-o",
            str(out),
            "-q",
            "--nomenclature",
            str(DERIVED),
            "--starts",
            "1",
        ],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr
    info = json.loads((out / "run.json").read_text())
    assert info["format"] == "ngsengine" and info["samples_used"] == 31
    assert info["loci_left_untyped"] == [{"id": "bc99", "reason": "A: result not approved (Rejected); left untyped"}]
    top = (out / "hfs.dat").read_text().splitlines()[0].split("\t")[0]
    assert top == "A*01:01~B*08:01"


def test_unknown_xml_layout_needs_format(tmp_path):
    p = tmp_path / "other.xml"
    p.write_text("<?xml version='1.0'?><Report><x/></Report>")
    run = subprocess.run([sys.executable, "-m", "hlahaplotools", "summary", str(p)], capture_output=True, text=True)
    assert run.returncode != 0 and "pass --format" in run.stderr


def test_missing_input_gives_a_clear_error(tmp_path):
    commands = (
        ["estimate", "--loci", "A:2f"],
        ["stats", "--loci", "A:2f"],
        ["families", "--loci", "A:2f"],
        ["summary"],
    )
    for cmd in commands:
        run = subprocess.run(
            [sys.executable, "-m", "hlahaplotools", *cmd, str(tmp_path / "nope.glc")], capture_output=True, text=True
        )
        assert run.returncode != 0, cmd
        assert "input not found" in run.stderr and "Traceback" not in run.stderr, (cmd, run.stderr)


def test_bad_locus_gives_one_line_error(tmp_path):
    p = tmp_path / "x.dat"
    p.write_text("id\tA\tA\tB\tB\ns1\t01:01\t02:01\t08:01\t07:02\n")
    run = subprocess.run(
        [sys.executable, "-m", "hlahaplotools", "summary", str(p), "--loci", "A,DPB1"], capture_output=True, text=True
    )
    assert run.returncode != 0
    assert run.stderr.startswith("error:") and "DPB1" in run.stderr and "Traceback" not in run.stderr, run.stderr
