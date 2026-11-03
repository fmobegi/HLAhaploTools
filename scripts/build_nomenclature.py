"""Download the reference data into this clone's data/ folder.

The same as `hlahaplotools download-data`, which works from any folder; this
wrapper keeps the clone's data/ as the default. Options:

    uv run python scripts/build_nomenclature.py                # pinned release
    uv run python scripts/build_nomenclature.py --release 3.65.0
    uv run python scripts/build_nomenclature.py --refresh-mac  # re-download MACs
    uv run python scripts/build_nomenclature.py --skip-mac     # IMGT/HLA files only
"""

from pathlib import Path

from hlahaplotools.reference_data import PINNED_RELEASE, main

__all__ = ["PINNED_RELEASE", "main"]

if __name__ == "__main__":
    main(default_data=Path(__file__).resolve().parent.parent / "data")
