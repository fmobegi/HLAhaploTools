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
    ds = h.read_mac(f, half_typed="unknown")
    assert ds.sample(0)["A"] == "A*01:01+"
    assert ds.sample(1)["A"] is None
    assert ds.sample(1)["B"] == "B*07:02+B*07:02"
    assert ds.missing_counts() == {"A": 1, "B": 0}
    assert ds.one_allele_counts() == {"A": 1, "B": 0}
    assert ds.fill_counts() == {}


def test_blank_second_allele_means_homozygous_by_default(tmp_path):
    f = tmp_path / "x.dat"
    f.write_text("id\tA\tA\tB\tB\n1\t01:01\t\t07:02\t08:01\n2\t\t\t07:02\t07:02\n")
    ds = h.read_mac(f)
    assert ds.sample(0)["A"] == "A*01:01+A*01:01"
    assert ds.sample(1)["A"] is None  # a fully blank locus stays untyped
    assert ds.one_allele_counts() == {"A": 0, "B": 0}
    assert ds.fill_counts()["A"]["homozygous"] == 1
    assert set(ds.fill_counts()) == {"A"}
    with pytest.raises(ValueError):
        h.read_mac(f, half_typed="maybe")


def test_drb345_follow_drb1(tmp_path):
    f = tmp_path / "x.dat"
    f.write_text(
        "id\tDRB1\tDRB1\tDRB3\tDRB3\tDRB4\tDRB4\tDRB5\tDRB5\n"
        # DR52 + DR53: one DRB3 and one DRB4 copy, no DRB5.
        "1\t03:01\t07:01\t01:01\t\t01:01\t\t\t\n"
        # Two DR52 groups: a single DRB3 allele is homozygous.
        "2\t03:01\t11:01\t02:02\t\t\t\t\t\n"
        # DRB1*01 and *08 carry none of the three genes.
        "3\t01:01\t08:01\t\t\t\t\t\t\n"
    )
    ds = h.read_mac(f)
    assert ds.sample(0)["DRB3"] == "DRB3*01:01+DRB3*NNNN"
    assert ds.sample(0)["DRB4"] == "DRB4*01:01+DRB4*NNNN"
    assert ds.sample(0)["DRB5"] == "DRB5*NNNN+DRB5*NNNN"
    assert ds.sample(1)["DRB3"] == "DRB3*02:02+DRB3*02:02"
    assert ds.sample(2)["DRB3"] == "DRB3*NNNN+DRB3*NNNN"
    fills = ds.fill_counts()
    assert fills["DRB3"] == {
        "homozygous": 1,
        "absent_gene": 2,
        "drb_completed": 1,
        "duplicate_read_as_one": 0,
        "drb1_conflicts": 0,
    }
    assert fills["DRB5"]["absent_gene"] == 6
    # A single DRB3/4/5 copy written twice is read as one copy.
    g = tmp_path / "y.dat"
    g.write_text(
        "id\tDRB1\tDRB1\tDRB3\tDRB3\tDRB4\tDRB4\n"
        "1\t01:01\t11:01\t02:02\t02:02\t\t\n"  # one DRB3 copy, written twice
        "2\t03:01\t11:01\t02:02\t02:02\t\t\n"  # two DRB3 copies: homozygous
        "3\t01:01\t08:01\t\t\t01:01\t01:01\n"  # DRB4 typed, DRB1 allows none
    )
    dup = h.read_mac(g)
    assert dup.sample(0)["DRB3"] == "DRB3*02:02+DRB3*NNNN"
    assert dup.sample(1)["DRB3"] == "DRB3*02:02+DRB3*02:02"
    assert dup.sample(2)["DRB4"] == "DRB4*01:01+DRB4*01:01"
    assert dup.fill_counts()["DRB3"]["duplicate_read_as_one"] == 1
    assert dup.fill_counts()["DRB4"]["drb1_conflicts"] == 1
    raw = h.read_mac(f, drb345=False)
    assert raw.sample(0)["DRB3"] == "DRB3*01:01+"
    assert raw.sample(2)["DRB3"] is None
