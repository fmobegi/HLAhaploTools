"""Read Hapl-o-Mat's example populations and the downloaded data.

Needs the reference clone and the downloaded data:
    git clone https://github.com/DKMS/Hapl-o-Mat reference/Hapl-o-Mat
    uv run python scripts/build_nomenclature.py
Tests that need them are skipped when they are missing.
"""

import json

import pytest
from conftest import EXAMPLES, MAC_TABLE, MANIFEST, needs_examples, needs_mac

import hlahaplotools as h

ALL_LOCI = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]


@needs_examples
@needs_mac
def test_mac_example_parses_fully():
    ds = h.read_mac(EXAMPLES / "populationData_a.dat", mac_codes=MAC_TABLE)
    assert ds.issues == []
    assert len(ds) == 100
    assert ds.loci == ALL_LOCI
    # Sample 2 has A*24:DXTW + A*02:570; DXTW expands to several A*24 alleles.
    left, right = ds.sample(1)["A"].split("+")
    assert len(left.split("/")) > 1 and all(a.startswith("A*24:") for a in left.split("/"))
    assert right == "A*02:570"


@needs_examples
def test_mac_example_without_table_reports_mac_codes():
    ds = h.read_mac(EXAMPLES / "populationData_a.dat")
    assert len(ds) + len(ds.issues) == 100
    assert all("MAC code" in reason for _, _, reason in ds.issues)


@needs_examples
def test_glsc_example_parses_fully():
    ds = h.read_glsc(EXAMPLES / "populationData_b.glc")
    assert ds.issues == []
    assert len(ds) == 100
    assert ds.loci == ALL_LOCI
    assert all(n == 0 for n in ds.missing_counts().values())


@pytest.mark.skipif(not MANIFEST.exists(), reason="run scripts/build_nomenclature.py first")
def test_manifest_records_release_mac_date_and_tables():
    m = json.loads(MANIFEST.read_text())
    assert m["imgthla"]["release"]
    assert set(m["imgthla"]["files"]) >= {"Allelelist.txt", "hla_nom_g.txt", "hla_nom_p.txt", "hla_nom.txt"}
    assert set(m["derived"]["rows"]) == {"AllAllelesExpanded.txt", "LargeG.txt", "P.txt", "Smallg.txt"}
    if "mac" in m:
        assert m["mac"]["last_updated"].startswith("LAST UPDATED")


def test_half_typed_and_missing_loci(tmp_path):
    f = tmp_path / "x.dat"
    f.write_text("id\tA\tA\tB\tB\n1\t01:01\t\t07:02\t08:01\n2\t\t\t07:02\t07:02\n")
    ds = h.read_mac(f)
    assert ds.sample(0)["A"] == "A*01:01+"
    assert ds.sample(1)["A"] is None
    assert ds.sample(1)["B"] == "B*07:02+B*07:02"
    assert ds.missing_counts() == {"A": 1, "B": 0}
    assert ds.one_allele_counts() == {"A": 1, "B": 0}
