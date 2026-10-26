"""Regression test on the real example pedigrees in examples/families.

The only pedigree problem in these files is Family3/C1, whose row is the
mother's typing (all loci but A) copied into the child's row.
"""

import sys

import pytest
from conftest import DERIVED, MAC_TABLE, ROOT, needs_derived, needs_mac

import hlahaplotools as h

EXAMPLES = ROOT / "examples" / "families"
pytestmark = [needs_derived, needs_mac]


@pytest.fixture(scope="module")
def cleaned(tmp_path_factory):
    """Both example tables, cleaned as validate_pedigrees.py does (leading zeros, uncovered loci)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from validate_pedigrees import clean, load_rows, write

    nom = h.Nomenclature(str(DERIVED))
    header, rows = load_rows(sorted(EXAMPLES.glob("*.csv")))
    loci, uncovered, _ = clean(header, rows, nom, "2f")
    path = tmp_path_factory.mktemp("examples") / "families.csv"
    write(path, header, rows)
    res = {l: ("asItIs" if l in uncovered else "2f") for l in loci}
    return h.read_family_table(path, mac_codes=MAC_TABLE), nom, res


def test_example_pedigrees_load(cleaned):
    data, _, _ = cleaned
    assert len(data) == 40
    assert len(data.families()) == 9
    assert data.issues == []


def test_only_family3_c1_is_flagged(cleaned):
    data, nom, res = cleaned
    flags = h.audit_families(data, nom, res)
    assert {f["sample"] for f in flags} == {"Family3/C1"}
    kinds = {f["kind"] for f in flags}
    assert {"parent_excluded", "possible_duplicate"} <= kinds
    dup = next(f for f in flags if f["kind"] == "possible_duplicate")
    assert dup["relative"] == "Family3/M"


def test_c1_defines_haplotypes_a_and_c_even_when_listed_later(nom):
    """In family_typing_data.csv, Family1 lists C2 first and Family2 lists C3 first.

    Haplotype labels must still be anchored on C1: C1 inherits A (father) and C (mother).
    """
    loci = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
    data = h.read_family_table(EXAMPLES / "family_typing_data.csv", loci=loci)
    phases = {p.family: p for p in h.phase_families(data, nom, "2f")}
    expected = {
        "Family1": [("Family1/C1", "A", "C"), ("Family1/C2", "B", "C"), ("Family1/C3", "A", "D")],
        "Family2": [("Family2/C1", "A", "C"), ("Family2/C2", "B", "D"), ("Family2/C3", "B", "D")],
    }
    for fam, children in expected.items():
        assert phases[fam].children == children
        labelled = phases[fam].labelled_haplotypes
        assert set(labelled) == {"A", "B", "C", "D"}
        assert labelled["A"] == phases[fam].haplotypes[0] and labelled["D"] == phases[fam].haplotypes[3]
    # Family1/C1 carries A*30:02 from the father and A*26:01 from the mother.
    assert phases["Family1"].labelled_haplotypes["A"].startswith("A*30:02~")
    assert phases["Family1"].labelled_haplotypes["C"].startswith("A*26:01~")


@pytest.fixture(scope="module")
def nom():
    return h.Nomenclature(str(DERIVED))


def test_examples_phase_without_inconsistencies_except_family3(cleaned):
    data, nom, res = cleaned
    phases = h.phase_families(data, nom, res)
    assert {p.family for p in phases if p.inconsistent_loci} == {"Family3"}
    assert all(p.status != "failed" for p in phases)
