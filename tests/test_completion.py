"""Filling in untyped loci (linked or exhaustive completion) and EM starts in the core."""

import pytest
from conftest import DERIVED, needs_derived

import hlahaplotools as h

pytestmark = needs_derived

# DQA1 travels with DQB1. People 5-8 are untyped at DQA1: linked completion should
# only try the DQA1 alleles that go with their DQB1 allele.
ROWS = [
    ("1", "01:01", "02:01", "05:01", "03:01", "02:01", "03:02"),
    ("2", "01:01", "01:01", "05:01", "05:01", "02:01", "02:01"),
    ("3", "02:01", "02:01", "03:01", "03:01", "03:02", "03:02"),
    ("4", "01:01", "02:01", "05:01", "03:01", "02:01", "03:02"),
    ("5", "03:01", "01:01", "", "", "02:01", "06:02"),
    ("6", "03:01", "03:01", "", "", "06:02", "06:02"),
    ("7", "02:01", "03:01", "", "", "03:02", "06:02"),
    ("8", "11:01", "03:01", "", "", "05:01", "06:02"),
    ("9", "03:01", "03:01", "01:02", "01:02", "06:02", "06:02"),
    ("10", "03:01", "11:01", "01:02", "01:01", "06:02", "05:01"),
]


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    p = tmp_path_factory.mktemp("c") / "x.dat"
    p.write_text("id\tA\tA\tDQA1\tDQA1\tDQB1\tDQB1\n" + "".join("\t".join(r) + "\n" for r in ROWS))
    return h.read_mac(str(p), half_typed="unknown"), h.Nomenclature(str(DERIVED))


def test_linked_tries_fewer_haplotypes_than_all(data):
    ds, nom = data
    linked = h.estimate(ds, nom, "2f", completion="linked", epsilon=1e-8)
    full = h.estimate(ds, nom, "2f", completion="all", epsilon=1e-8)
    assert linked.n_samples == full.n_samples == len(ROWS)
    assert linked.n_completions <= full.n_completions
    assert abs(sum(f for _, f in linked.haplotypes) - 1) < 1e-9
    # The common haplotype is estimated alike either way.
    top = dict(full.haplotypes)
    for hap, f in linked.haplotypes[:2]:
        assert abs(f - top.get(hap, 0)) < 0.05
    with pytest.raises(ValueError):
        h.estimate(ds, nom, "2f", completion="some")


def test_starts_run_in_the_core(data):
    ds, nom = data
    est = h.estimate(ds, nom, "2f", starts=4, seed=3, epsilon=1e-8)
    assert [(i, s) for i, s, *_ in est.starts] == [
        ("numberOccurrence", 3),
        ("perturbation", 4),
        ("perturbation", 5),
        ("perturbation", 6),
    ]
    best = max(ll for _, _, ll, _, _ in est.starts)
    assert est.log_likelihood >= best - 1e-9 * abs(best)
    one = h.estimate(ds, nom, "2f", epsilon=1e-8)
    assert len(one.starts) == 1
    with pytest.raises(ValueError):
        h.estimate(ds, nom, "2f", starts=0)


def test_diplotypes_can_be_skipped(data):
    ds, nom = data
    assert len(h.estimate(ds, nom, "2f", diplotypes=False).best_diplotypes()) == 0
    assert len(h.estimate(ds, nom, "2f").best_diplotypes()) == len(ROWS)
