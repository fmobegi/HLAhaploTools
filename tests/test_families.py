"""Family data: reading, Mendel check, phasing, founder-only estimation, comparison."""

import json
import subprocess
import sys

import pytest
from conftest import DERIVED, ROOT, needs_derived

import hlahaplotools as h

pytestmark = needs_derived


@pytest.fixture(scope="module")
def nom():
    return h.Nomenclature(str(DERIVED))


@pytest.fixture(scope="module")
def sim(tmp_path_factory):
    """200 simulated families (1-4 children, 10% of parents untyped) with known founders."""
    out = tmp_path_factory.mktemp("families")
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "simulate.py"),
            str(out),
            "--format",
            "families",
            "--samples",
            "200",
            "--children",
            "1-4",
            "--missing-parent",
            "0.1",
            "--seed",
            "5",
        ],
        check=True,
        capture_output=True,
    )
    truth = {}
    for line in (out / "truth_founders.tsv").read_text().splitlines():
        fam, role, a, b = line.split("\t")
        truth[(fam, role)] = sorted([a, b])
    return out, truth


def test_reads_families_and_roles(sim):
    out, _ = sim
    data = h.read_family_table(out / "families.csv")
    fams = data.families()
    assert len(fams) == 200
    assert data.unrelated == [] and data.issues == []
    assert all(f["children"] for f in fams)
    assert data.loci == ["A", "B", "C", "DRB1", "DQB1"]


def test_phasing_matches_truth(sim, nom):
    out, truth = sim
    phases = h.phase_families(h.read_family_table(out / "families.csv"), nom, "asItIs")
    checked = 0
    for p in phases:
        for i, role in enumerate("FM"):
            pair = p.haplotypes[2 * i : 2 * i + 2]
            if p.status != "failed" and not any("*?" in x for x in pair):
                assert sorted(pair) == truth[(p.family, role)], (p.family, role, p.notes)
                checked += 1
    assert checked > 300  # most parents are fully phased
    assert sum(p.status == "complete" for p in phases) > 140


def test_mendel_error_is_found(sim, nom, tmp_path):
    out, _ = sim
    lines = (out / "families.csv").read_text().splitlines()
    # Give the first child an A allele neither parent has.
    i = next(i for i, l in enumerate(lines) if ",C1," in l)
    cells = lines[i].split(",")
    cells[2] = "A*99:99"
    lines[i] = ",".join(cells)
    bad = tmp_path / "bad.csv"
    bad.write_text("\n".join(lines) + "\n")
    issues = h.check_mendelian(h.read_family_table(bad), nom, "asItIs")
    assert any(m["locus"] == "A" and m["sample"] == f"{cells[0]}/C1" for m in issues)


def test_founder_estimate_and_comparison(sim, nom):
    out, truth = sim
    data = h.read_family_table(out / "families.csv")
    est, phases = h.estimate_families(data, nom, "asItIs", epsilon=1e-8)
    assert est.converged
    assert 380 <= est.n_samples <= 400  # two founders per family, typed or rebuilt
    assert abs(sum(f for _, f in est.haplotypes) - 1.0) < 1e-9
    dips = {d[0]: d for d in est.best_diplotypes()}
    complete = [p for p in phases if p.status == "complete"]
    for p in complete[:20]:
        sid, h1, h2, prob = dips[p.parent_id(0)]
        assert sorted([h1, h2]) == sorted(p.haplotypes[:2]) and prob == pytest.approx(1.0)

    unphased, _ = h.estimate_families(data, nom, "asItIs", use_family_phase=False, epsilon=1e-8)
    cmp = h.compare_to_families(unphased, phases)
    assert cmp["n_compared"] > 250
    assert cmp["n_exact"] / cmp["n_compared"] > 0.9
    assert 0.0 <= cmp["total_variation"] < 0.2
    # Phased input reproduces the family haplotypes exactly.
    assert h.compare_to_families(est, phases)["n_exact"] == h.compare_to_families(est, phases)["n_compared"]


def test_audit_finds_planted_errors(nom, tmp_path):
    """Planted non-paternity, copied rows, recombinants and typos are flagged; clean children rarely are."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "simulate.py"),
            str(tmp_path),
            "--format",
            "families",
            "--samples",
            "600",
            "--children",
            "2-4",
            "--loci",
            "A,C,B,DRB1,DQB1,DPB1",
            "--haplotypes",
            "1500",
            "--alleles",
            "40",
            "--seed",
            "9",
            "--errors",
            "other_father:0.05,duplicate_of_mother:0.03,recombinant:0.05,typo:0.05",
        ],
        check=True,
        capture_output=True,
    )
    planted = dict(line.split("\t")[1:] for line in (tmp_path / "truth_errors.tsv").read_text().splitlines())
    data = h.read_family_table(tmp_path / "families.csv")
    flags = h.audit_families(data, nom, "asItIs")
    flagged = {}
    for f in flags:
        if f["sample"]:
            flagged.setdefault(f["sample"], set()).add(f["kind"])
    children = [c for fam in data.families() for c in fam["children"]]
    expected = {
        "other_father": "parent_excluded",
        "duplicate_of_mother": "possible_duplicate",
        "recombinant": "recombination",
        "typo": "single_locus_mismatch",
    }
    for kind, flag in expected.items():
        hits = [c for c, k in planted.items() if k == kind]
        found = sum(flag in flagged.get(c, set()) for c in hits)
        assert hits and found / len(hits) > 0.6, (kind, found, len(hits))
    clean = [c for c in children if c not in planted]
    false_alarms = sum(bool(flagged.get(c)) for c in clean)
    assert false_alarms / len(clean) < 0.05, false_alarms


def test_families_cli(sim, tmp_path):
    out, _ = sim
    res = tmp_path / "res"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "hlahaplotools",
            "families",
            str(out / "families.csv"),
            "--loci",
            "A:asItIs,B:asItIs,C:asItIs,DRB1:asItIs,DQB1:asItIs",
            "-o",
            str(res),
            "-q",
            "--data-dir",
            str(ROOT / "data"),
        ],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr
    for name in [
        "flags.tsv",
        "mendel.tsv",
        "phases.tsv",
        "children.tsv",
        "hfs.dat",
        "diplotypes.tsv",
        "comparison.json",
    ]:
        assert (res / name).exists(), name
    info = json.loads((res / "run.json").read_text())
    assert info["families"] == 200
    # No per-locus Mendel failures in clean data; a few families with no typed
    # parents and a single child are reported as unphasable, which is correct.
    mendel_rows = [l.split("\t") for l in (res / "mendel.tsv").read_text().splitlines()[1:]]
    assert all(row[2] == "" for row in mendel_rows)
    assert info["phasing"]["failed"] == len(mendel_rows) <= 5
