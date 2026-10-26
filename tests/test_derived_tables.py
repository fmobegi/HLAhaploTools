"""The lookup tables built by scripts/derive_tables.py must equal Hapl-o-Mat's own."""

import pytest
from conftest import HAPLOMAT_DIR, IMGT_RELEASE_DIR

TABLES = ["AllAllelesExpanded.txt", "LargeG.txt", "P.txt", "Smallg.txt"]


@pytest.mark.skipif(
    IMGT_RELEASE_DIR is None or not (HAPLOMAT_DIR / "data" / "LargeG.txt").exists(),
    reason="needs data/imgthla and Hapl-o-Mat's data/ built from the same release",
)
def test_tables_identical_to_haplomat_prepare_data(tmp_path):
    from derive_tables import derive_all

    derive_all(IMGT_RELEASE_DIR, tmp_path)
    for name in TABLES:
        ours = (tmp_path / name).read_bytes()
        theirs = (HAPLOMAT_DIR / "data" / name).read_bytes()
        assert ours == theirs, f"{name} differs from Hapl-o-Mat's prepareData output"
