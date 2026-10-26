"""Shared paths and skip markers for tests that need downloaded data."""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))  # for derive_tables

HAPLOMAT_DIR = ROOT / "reference" / "Hapl-o-Mat"
EXAMPLES = HAPLOMAT_DIR / "examplePopulations"
MANIFEST = ROOT / "data" / "manifest.json"
MAC_TABLE = ROOT / "data" / "mac" / "alpha.v3.txt"

RELEASE = json.loads(MANIFEST.read_text())["imgthla"]["release"] if MANIFEST.exists() else None
IMGT_RELEASE_DIR = ROOT / "data" / "imgthla" / RELEASE if RELEASE else None
DERIVED = ROOT / "data" / "derived" / RELEASE if RELEASE else None

needs_examples = pytest.mark.skipif(not EXAMPLES.exists(), reason="Hapl-o-Mat clone not in reference/")
needs_mac = pytest.mark.skipif(not MAC_TABLE.exists(), reason="run scripts/build_nomenclature.py first")
needs_derived = pytest.mark.skipif(
    DERIVED is None or not DERIVED.exists(), reason="run scripts/build_nomenclature.py first"
)
needs_haplomat = pytest.mark.skipif(
    not (HAPLOMAT_DIR / "haplomat").exists() or not (HAPLOMAT_DIR / "data" / "LargeG.txt").exists(),
    reason="build reference/Hapl-o-Mat (make) and its data/ first",
)
