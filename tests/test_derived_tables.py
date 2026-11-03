"""The lookup tables built by scripts/derive_tables.py must equal Hapl-o-Mat's own."""

import pytest
from conftest import HAPLOMAT_DIR, IMGT_RELEASE_DIR

TABLES = ["AllAllelesExpanded.txt", "LargeG.txt", "P.txt", "Smallg.txt"]


@pytest.mark.skipif(
    IMGT_RELEASE_DIR is None or not (HAPLOMAT_DIR / "data" / "LargeG.txt").exists(),
    reason="needs data/imgthla and Hapl-o-Mat's data/ built from the same release",
)
def test_tables_identical_to_haplomat_prepare_data(tmp_path):
    from hlahaplotools.derive_tables import derive_all

    derive_all(IMGT_RELEASE_DIR, tmp_path, extra_loci=False)
    for name in TABLES:
        ours = (tmp_path / name).read_bytes()
        theirs = (HAPLOMAT_DIR / "data" / name).read_bytes()
        assert ours == theirs, f"{name} differs from Hapl-o-Mat's prepareData output"


@pytest.mark.skipif(
    IMGT_RELEASE_DIR is None or not (HAPLOMAT_DIR / "data" / "LargeG.txt").exists(),
    reason="needs data/imgthla and Hapl-o-Mat's data/ built from the same release",
)
def test_extra_loci_only_add_rows(tmp_path):
    """MICA, MICB, H, J, ... are added; every row for Hapl-o-Mat's loci is unchanged."""
    from hlahaplotools.derive_tables import derive_all

    derive_all(IMGT_RELEASE_DIR, tmp_path, extra_loci=True)
    for name in TABLES:
        ours = (tmp_path / name).read_text().splitlines()
        theirs = (HAPLOMAT_DIR / "data" / name).read_text().splitlines()
        theirs_loci = {l.split("*", 1)[0] for l in theirs}
        assert [l for l in ours if l.split("*", 1)[0] in theirs_loci] == theirs, name
    expanded = {
        l.split("	", 1)[0]: l.split("	")[1:]
        for l in (tmp_path / "AllAllelesExpanded.txt").read_text().splitlines()
    }
    for locus in ("MICA", "MICB", "H", "J"):
        assert any(k.startswith(locus + "*") for k in expanded), locus
    assert "MICA*008:01:01" in expanded["MICA*008:01"]
