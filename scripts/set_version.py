"""Set the version everywhere it is written, so the files never disagree.

    python scripts/set_version.py 0.1.0rc3      # a release candidate
    python scripts/set_version.py 0.1.0         # a release

Updates pyproject.toml (the package; PEP 440, e.g. 0.1.0rc3), both Cargo.toml files (the Rust crates;
SemVer, e.g. 0.1.0-rc.3), and the version in .zenodo.json and CITATION.cff (e.g. 0.1.0-rc3, matching the
release tag v0.1.0-rc3). It does not touch DOIs or dates: put the new DOI in CITATION.cff and the README
once Zenodo has minted it. Cargo.lock and uv.lock follow on the next `cargo build` or `uv sync`.
The TestPyPI workflow also calls this for its .devN builds.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PEP440 = re.compile(r"(\d+\.\d+\.\d+)(?:(a|b|rc)(\d+))?(?:\.dev(\d+))?")


def forms(version: str) -> tuple[str, str]:
    """(SemVer for Cargo, tag-style for Zenodo and CITATION.cff) from a PEP 440 version."""
    m = PEP440.fullmatch(version)
    if not m:
        sys.exit(f"error: {version!r} is not a version like 0.1.0, 0.1.0rc3 or 0.1.0.dev5")
    base, pre, num, dev = m.groups()
    semver = base + (f"-{pre}.{num}" if pre else "") + (f"{'.' if pre else '-'}dev.{dev}" if dev else "")
    tag = base + (f"-{pre}{num}" if pre else "") + (f".dev{dev}" if dev else "")
    return semver, tag


def replace(path: Path, pattern: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    out, n = re.subn(pattern, new, text, count=1, flags=re.M)
    if n != 1:
        sys.exit(f"error: no version found in {path.relative_to(ROOT)}")
    path.write_text(out, encoding="utf-8", newline="\n")


def main() -> None:
    version = sys.argv[1]
    semver, tag = forms(version)
    replace(ROOT / "pyproject.toml", r'^version = "[^"]*"$', f'version = "{version}"')
    for cargo in (ROOT / "Cargo.toml", ROOT / "crates" / "hlahaplotools-core" / "Cargo.toml"):
        replace(cargo, r'^version = "[^"]*"$', f'version = "{semver}"')
    replace(ROOT / ".zenodo.json", r'"version": "[^"]*"', f'"version": "{tag}"')
    replace(ROOT / "CITATION.cff", r'^version: "[^"]*"$', f'version: "{tag}"')
    print(f"version {version} (Cargo {semver}, release {tag})")


if __name__ == "__main__":
    main()
