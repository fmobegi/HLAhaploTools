"""EM from several start values."""

import json
import subprocess
import sys

import pytest
from conftest import DERIVED, ROOT, needs_derived

import hlahaplotools as h
from hlahaplotools.multistart import estimate_multistart

POPULATION = ROOT / "examples" / "population" / "population.dat"
LOCI = {"A": "2f", "C": "2f", "B": "2f", "DRB1": "2f", "DQB1": "2f"}
pytestmark = needs_derived


@pytest.fixture(scope="module")
def data():
    return h.read_mac(str(POPULATION)), h.Nomenclature(str(DERIVED))


def test_best_of_several_starts(data):
    ds, nom = data
    best, summary = estimate_multistart(ds, nom, LOCI, starts=4, seed=7, epsilon=1e-8)
    assert [(i, s) for i, s, *_ in summary] == [
        ("numberOccurrence", 7),
        ("perturbation", 8),
        ("perturbation", 9),
        ("perturbation", 10),
    ]
    top = max(ll for _, _, ll, _, _ in summary)
    assert best.log_likelihood >= top - 1e-9 * abs(top)  # near-ties keep the earlier start
    # Reproducible: the same call gives the same answer.
    again, _ = estimate_multistart(ds, nom, LOCI, starts=4, seed=7, epsilon=1e-8)
    assert again.haplotypes == best.haplotypes


def test_one_start_is_the_plain_estimate(data):
    ds, nom = data
    best, summary = estimate_multistart(ds, nom, LOCI, starts=1, epsilon=1e-8)
    plain = h.estimate(ds, nom, LOCI, epsilon=1e-8)
    assert len(summary) == 1 and best.haplotypes == plain.haplotypes
    with pytest.raises(ValueError):
        estimate_multistart(ds, nom, LOCI, starts=0)


def test_cli_starts_are_recorded(tmp_path):
    loci = ",".join(f"{l}:{r}" for l, r in LOCI.items())
    run = subprocess.run(
        [sys.executable, "-m", "hlahaplotools", "estimate", str(POPULATION), "--loci", loci, "-o", str(tmp_path),
         "--starts", "3", "--data-dir", str(ROOT / "data")],
        capture_output=True, text=True,
    )  # fmt: skip
    assert run.returncode == 0, run.stderr
    assert "Starts: 3; " in run.stdout
    info = json.loads((tmp_path / "run.json").read_text())
    assert info["parameters"]["starts"] == 3
    assert len(info["starts"]) == 3
    top = max(s["log_likelihood"] for s in info["starts"])
    assert info["log_likelihood"] >= top - 1e-9 * abs(top)
