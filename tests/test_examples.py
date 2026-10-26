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


def test_examples_phase_without_inconsistencies_except_family3(cleaned):
    data, nom, res = cleaned
    phases = h.phase_families(data, nom, res)
    assert {p.family for p in phases if p.inconsistent_loci} == {"Family3"}
    assert all(p.status != "failed" for p in phases)
