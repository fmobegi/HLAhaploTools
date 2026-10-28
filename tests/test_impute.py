"""Imputation from known frequencies, the genotype list, serology, V2 names, phased GL strings."""

import json
import subprocess
import sys

import pytest
from conftest import DERIVED, ROOT, needs_derived

import hlahaplotools as h

POP = ROOT / "examples" / "population" / "population.dat"


def run(*args):
    r = subprocess.run([sys.executable, "-m", "hlahaplotools", *map(str, args)], capture_output=True, text=True)
    assert r.returncode in (0, 2), r.stderr
    return r


@needs_derived
def test_impute_ranks_pairs_that_fit_and_sums_to_one(tmp_path):
    nom = h.Nomenclature(str(DERIVED))
    ds = h.read_mac(str(POP), loci=["A", "B"])
    est = h.estimate(ds, nom, "g", epsilon=1e-8)
    rows, issues = h.impute(ds, nom, [("all", est.haplotypes, 1.0)], "g", top=1000)
    assert len(rows) + len(issues) == len(ds) - len(est.issues)
    for r in rows[:50]:
        assert abs(sum(p[3] for p in r["pairs"]) - 1) < 1e-9 or r["n_pairs"] > 1000
        assert abs(sum(g[1] for g in r["genotypes"]) - 1) < 1e-9 or len(r["genotypes"]) == 1000
        assert r["pairs"] == sorted(r["pairs"], key=lambda p: -p[3])
    # The best pair agrees with the estimate's own most likely diplotype.
    best = {sid: {a, b} for sid, a, b, _ in est.best_diplotypes()}
    agree = sum({r["pairs"][0][1], r["pairs"][0][2]} == best[r["id"]] for r in rows)
    assert agree / len(rows) > 0.99


@needs_derived
def test_impute_command_takes_loci_from_run_json(tmp_path):
    est_dir = tmp_path / "est"
    run("estimate", POP, "--loci", "A:g,C:g,B:g", "-o", est_dir, "--nomenclature", DERIVED, "-q", "--starts", "1")
    out = tmp_path / "imp"
    run("impute", POP, "--frequencies", est_dir / "hfs.dat", "-o", out, "--nomenclature", DERIVED, "-q", "--top", "3")
    pairs = (out / "pairs.tsv").read_text().splitlines()
    assert pairs[0].split("\t") == ["id", "rank", "population", "haplotype_1", "haplotype_2", "probability"]
    assert pairs[1].split("\t")[3].count("~") == 2
    info = json.loads((out / "run.json").read_text())
    assert info["loci_and_resolutions"] == {"A": "g", "C": "g", "B": "g"} and info["samples_imputed"] > 0
    # Two populations: posteriors are written per person.
    out2 = tmp_path / "imp2"
    run("impute", POP, "--frequencies", f"x={est_dir / 'hfs.dat'}", "--frequencies", f"y={est_dir / 'hfs.dat'}",
        "--prior", "x=3", "-o", out2, "--nomenclature", DERIVED, "-q")  # fmt: skip
    post = [l.split("\t") for l in (out2 / "populations.tsv").read_text().splitlines()[1:3]]
    assert [p[1] for p in post] == ["x", "y"] and abs(float(post[0][2]) - 0.75) < 1e-9


@needs_derived
def test_genotype_list_round_trip(tmp_path):
    """--write-genotypes, then --format genotypes with asItIs, gives the same estimate."""
    a = tmp_path / "a"
    run("estimate", POP, "--loci", "A:g,B:g", "-o", a, "--nomenclature", DERIVED, "-q", "--starts", "1",
        "--write-genotypes", "--epsilon", "1e-8")  # fmt: skip
    lines = (a / "genotypes.dat").read_text().splitlines()
    sid, kind, w, code = lines[0].split("\t")
    assert len(kind) == 2 and 0 < float(w) <= 1 and code.count("^") == 1
    b = tmp_path / "b"
    run("estimate", a / "genotypes.dat", "--format", "genotypes", "--loci", "A:asItIs,B:asItIs", "-o", b,
        "--nomenclature", DERIVED, "-q", "--starts", "1", "--epsilon", "1e-8")  # fmt: skip
    fa = dict(l.split("\t") for l in (a / "hfs.dat").read_text().splitlines())
    fb = dict(l.split("\t") for l in (b / "hfs.dat").read_text().splitlines())
    assert set(fa) == set(fb)
    assert max(abs(float(fa[k]) - float(fb[k])) for k in fa) < 1e-6


@needs_derived
def test_serology_and_v2_names_are_explained_or_read(tmp_path):
    p = tmp_path / "s.dat"
    p.write_text("id\tA\tA\tB\tB\n1\t0201\t010101\t62\t08:01\n")
    ds = h.read_mac(str(p))
    assert not ds.issues  # 010101 is kept for the resolver, as a possible pre-2010 name
    plain = h.Nomenclature(str(DERIVED))
    with pytest.raises(ValueError, match="serological antigen B62"):
        plain.translate("B*62", "2f")
    with pytest.raises(ValueError, match="pre-2010 allele name .now A.01:01:01"):
        plain.translate("A*010101", "2f")
    with pytest.raises(ValueError, match="no usable samples"):
        h.estimate(ds, plain, "2f")
    full = h.Nomenclature(str(DERIVED), serology=True, v2_names=True)
    est = h.estimate(ds, full, "2f")
    assert est.n_samples == 1
    names = {x for hap, _ in est.haplotypes for x in hap.split("~")}
    assert {"A*02:01", "A*01:01"} <= names and any(n.startswith("B*15:") for n in names)
    with pytest.raises(ValueError, match="A*92"):
        plain.translate("A*92", "g")


@needs_derived
def test_phased_gl_string_is_read_without_phase(tmp_path):
    p = tmp_path / "p.glc"
    p.write_text("1\tHLA-A*01:01~HLA-B*08:01+HLA-A*02:01~HLA-B*07:02\n")
    ds = h.read_glsc(str(p))
    assert ds.loci == ["A", "B"] and len(ds) == 1
    assert any("without its phase" in n[2] for n in ds.notes)


def test_lgx_resolution_and_aliases():
    if DERIVED is None or not DERIVED.exists():
        pytest.skip("run scripts/build_nomenclature.py first")
    nom = h.Nomenclature(str(DERIVED))
    assert nom.translate("A*01:01:01:01", "lg") == nom.translate("A*01:01:01:01", "g")
    assert nom.translate("A*01:01:01:01", "lgx") == [("A*01:01", 1.0)]
