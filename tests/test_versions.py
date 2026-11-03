"""Every file that names the version must agree (set them together with scripts/set_version.py)."""

import json
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from set_version import forms  # noqa: E402


def test_versions_agree():
    version = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    semver, tag = forms(version)
    for cargo in (ROOT / "Cargo.toml", ROOT / "crates" / "hlahaplotools-core" / "Cargo.toml"):
        assert tomllib.loads(cargo.read_text(encoding="utf-8"))["package"]["version"] == semver, cargo
    assert json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))["version"] == tag
    cff = re.search(r'^version: "([^"]*)"$', (ROOT / "CITATION.cff").read_text(encoding="utf-8"), re.M)
    assert cff and cff.group(1) == tag


def test_version_forms():
    assert forms("0.1.0") == ("0.1.0", "0.1.0")
    assert forms("0.1.0rc2") == ("0.1.0-rc.2", "0.1.0-rc2")
    assert forms("0.1.0rc2.dev3") == ("0.1.0-rc.2.dev.3", "0.1.0-rc2.dev3")
    assert forms("0.1.0.dev4") == ("0.1.0-dev.4", "0.1.0.dev4")
