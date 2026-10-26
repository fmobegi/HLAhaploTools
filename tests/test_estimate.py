"""Estimation: API behaviour, missing data, determinism, CLI, and parity with Hapl-o-Mat."""

import json
import subprocess
import sys

import pytest
from conftest import DERIVED, EXAMPLES, ROOT, needs_derived, needs_examples, needs_haplomat

import hlahaplotools as h


@pytest.fixture(scope="module")
def nom():
    return h.Nomenclature(str(DERIVED))


def write(tmp_path, text, name="pop.dat"):
    p = tmp_path / name
    p.write_text(text)
    return p


def freq(est, hap):
    return dict(est.haplotypes).get(hap, 0.0)


@needs_derived
def test_translate_known_groups(nom):
    # A*01:01 at g resolution maps to its g group with full weight.
    (code, w), *rest = nom.translate("A*01:01", "g")
    assert code.startswith("A*01:01") and w > 0.5
    assert abs(sum(w for _, w in nom.translate("A*02", "2f")) - 1.0) < 1e-12
    with pytest.raises(ValueError):
        nom.translate("A*99:99", "g")
    with pytest.raises(ValueError):
        nom.translate("A*01:01", "6d")


@needs_derived
def test_estimate_simple_linkage(nom, tmp_path):
    rows = ["id\tA\tA\tB\tB"]
    rows += [f"h{i}\t01:01\t01:01\t08:01\t08:01" for i in range(30)]
    rows += [f"k{i}\t02:01\t02:01\t07:02\t07:02" for i in range(10)]
    rows += ["x\t01:01\t02:01\t08:01\t07:02"]
    est = h.estimate(h.read_mac(write(tmp_path, "\n".join(rows) + "\n")), nom, "2f", epsilon=1e-9)
    assert est.converged
    top = est.frequencies(cut=1e-3)
    assert top[0][0] == "A*01:01~B*08:01"
    assert abs(sum(f for _, f in est.frequencies(cut=0.0, renormalize=False)) - 1.0) < 1e-9


@needs_derived
def test_resolutions_as_dict_list_or_string(nom, tmp_path):
    ds = h.read_mac(write(tmp_path, "id\tA\tA\tB\tB\n1\t01:01\t02:01\t08:01\t07:02\n"))
    a = h.estimate(ds, nom, {"A": "g", "B": "g"})
    b = h.estimate(ds, nom, ["g", "g"])
    c = h.estimate(ds, nom, "g")
    d = h.estimate(ds, nom)
    assert a.haplotypes == b.haplotypes == c.haplotypes == d.haplotypes
    with pytest.raises(ValueError):
        h.estimate(ds, nom, {"A": "g"})


@needs_derived
def test_missing_loci_are_used_not_dropped(nom, tmp_path):
    rows = ["id\tA\tA\tB\tB"]
    rows += [f"h{i}\t01:01\t01:01\t08:01\t08:01" for i in range(10)]
    rows += [f"k{i}\t02:01\t02:01\t07:02\t07:02" for i in range(10)]
    rows += ["m1\t01:01\t01:01\t\t", "m2\t\t\t07:02\t07:02", "m3\t01:01\t\t08:01\t08:01"]
    est = h.estimate(h.read_mac(write(tmp_path, "\n".join(rows) + "\n")), nom, "2f", epsilon=1e-10)
    assert est.n_samples == 23 and est.n_samples_incomplete == 3
    # m1 adds two copies of A*01:01~B*08:01, m2 two of A*02:01~B*07:02,
    # m3 at least one more A*01:01~B*08:01.
    assert freq(est, "A*01:01~B*08:01") > 22 / 46
    assert abs(sum(f for _, f in est.haplotypes) - 1.0) < 1e-9


@needs_derived
def test_unknown_allele_leaves_sample_out(nom, tmp_path):
    ds = h.read_mac(write(tmp_path, "id\tA\tA\n1\t01:01\t02:01\n2\t99:99\t01:01\n"))
    est = h.estimate(ds, nom, "g")
    assert est.n_samples == 1
    assert est.issues[0][1] == "2" and "A*99:99" in est.issues[0][2]


@needs_derived
def test_threads_do_not_change_results(nom, tmp_path):
    rows = ["id\tA\tA\tB\tB\tC\tC"]
    a, b, c = ["01:01", "02:01", "03:01", "11:01"], ["07:02", "08:01", "44:02", "35:01"], ["07:01", "04:01", "05:01"]
    for i in range(2000):
        cc = ("", "") if i % 7 == 0 else (c[i % 3], c[(i // 3) % 3])
        rows.append(f"s{i}\t{a[i % 4]}\t{a[(i // 4) % 4]}\t{b[(i * 3) % 4]}\t{b[(i // 5) % 4]}\t{cc[0]}\t{cc[1]}")
    ds = h.read_mac(write(tmp_path, "\n".join(rows) + "\n"))
    one = h.estimate(ds, nom, "g", threads=1, epsilon=1e-8)
    many = h.estimate(ds, nom, "g", threads=4, epsilon=1e-8)
    assert one.haplotypes == many.haplotypes
    assert one.history == many.history


@needs_derived
def test_epsilon_larger_than_half_over_n_is_refused(nom, tmp_path):
    ds = h.read_mac(write(tmp_path, "id\tA\tA\n1\t01:01\t02:01\n2\t01:01\t03:01\n"))
    with pytest.raises(ValueError, match="epsilon"):
        h.estimate(ds, nom, "g", epsilon=0.3)


@needs_derived
@needs_examples
def test_cli_estimate_writes_outputs(tmp_path):
    out = tmp_path / "res"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "hlahaplotools",
            "estimate",
            str(EXAMPLES / "populationData_b.glc"),
            "--loci",
            "A:g,B:g,DRB1:g",
            "-o",
            str(out),
            "-q",
            "--data-dir",
            str(ROOT / "data"),
        ],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr
    lines = (out / "hfs.dat").read_text().splitlines()
    assert lines and all(len(l.split("\t")) == 2 for l in lines)
    info = json.loads((out / "run.json").read_text())
    assert info["imgthla_release"] and info["samples_used"] == 100 and info["converged"]
    assert (out / "epsilon.dat").exists()


@needs_haplomat
@needs_examples
@pytest.mark.parametrize(
    "fmt, data, loci",
    [
        ("MAC", "populationData_a.dat", "A:g,B:g,DRB1:1f"),
        ("GLSC", "populationData_b.glc", "A:g,C:g,DPB1:g,DQB1:g,DRB1:P"),
    ],
)
def test_matches_original_haplomat(fmt, data, loci):
    run = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "compare_haplomat.py"),
            fmt,
            str(EXAMPLES / data),
            "--loci",
            loci,
            "--tolerance",
            "1e-10",
        ],
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    assert "largest difference" in run.stdout
