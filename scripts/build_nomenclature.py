"""Download a pinned IMGT/HLA release and the NMDP MAC table into data/.

Usage (from the project root):

    uv run python scripts/build_nomenclature.py                # pinned release
    uv run python scripts/build_nomenclature.py --release 3.65.0
    uv run python scripts/build_nomenclature.py --refresh-mac  # re-download MACs

Layout written:

    data/
      imgthla/<release>/Allelelist.txt, hla_nom.txt, hla_nom_g.txt,
                        hla_nom_p.txt, Deleted_alleles.txt
      derived/<release>/AllAllelesExpanded.txt, LargeG.txt, P.txt, Smallg.txt
      mac/alpha.v3.txt
      manifest.json      # release, source URLs, sha256, MAC "LAST UPDATED"

The derived tables are what the resolver reads; see scripts/derive_tables.py.

IMGT/HLA files come from the ANHIG/IMGTHLA GitHub tag for the release, and
each file's "# version:" header is checked against it. The MAC table has no
release number, so it is downloaded once and kept; the manifest records its
"LAST UPDATED" line and checksum so results can be traced back to it.

Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import urllib.request
import zipfile
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_tables import derive_all  # noqa: E402

# Bump deliberately: every result file should name the release it used.
PINNED_RELEASE = "3.66.0"

IMGT_RAW = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/v{release}-alpha/{path}"
IMGT_FILES = {
    "Allelelist.txt": "Allelelist.txt",
    "hla_nom.txt": "wmda/hla_nom.txt",  # includes deleted and renamed alleles
    "hla_nom_g.txt": "wmda/hla_nom_g.txt",
    "hla_nom_p.txt": "wmda/hla_nom_p.txt",
    "Deleted_alleles.txt": "Deleted_alleles.txt",
}
MAC_URL = "https://hml.nmdp.org/mac/files/alpha.v3.zip"

ROOT = Path(__file__).resolve().parent.parent


def rel(path: Path) -> str:
    """Path relative to the project root when possible, for messages and the manifest."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    """Download to a temporary file first, so an interrupted run leaves nothing half-written."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "hlahaplotools-build"})
    print(f"  downloading {url}")
    with urllib.request.urlopen(req, timeout=120) as resp, tmp.open("wb") as out:
        shutil.copyfileobj(resp, out, length=1 << 20)
    tmp.replace(dest)
    print(f"    -> {rel(dest)} ({dest.stat().st_size / 1e6:.1f} MB)")


def imgt_version(path: Path) -> str | None:
    """Read the '# version: IPD-IMGT/HLA X.Y.Z' header line."""
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("#"):
                break
            m = re.match(r"#\s*version:\s*IPD-IMGT/HLA\s+(\S+)", line)
            if m:
                return m.group(1)
    return None


def fetch_imgt(release: str, data: Path, force: bool) -> dict:
    out_dir = data / "imgthla" / release
    files = {}
    for name, repo_path in IMGT_FILES.items():
        dest = out_dir / name
        url = IMGT_RAW.format(release=release, path=repo_path)
        if force or not dest.exists():
            download(url, dest)
        else:
            print(f"  have {rel(dest)}")
        found = imgt_version(dest)
        if found != release:
            sys.exit(f"error: {dest} says version {found!r}, expected {release!r}")
        files[name] = {"url": url, "sha256": sha256(dest), "bytes": dest.stat().st_size}
    return {"release": release, "dir": rel(out_dir), "files": files}


def fetch_mac(data: Path, refresh: bool) -> dict:
    mac_dir = data / "mac"
    txt = mac_dir / "alpha.v3.txt"
    if refresh or not txt.exists():
        zpath = mac_dir / "alpha.v3.zip"
        download(MAC_URL, zpath)
        with zipfile.ZipFile(zpath) as z:
            member = next(n for n in z.namelist() if n.endswith("alpha.v3.txt"))
            with z.open(member) as src, txt.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=1 << 20)
        zpath.unlink()  # the zip is ~110 MB; the text file is what we use
        print(f"    -> {rel(txt)} ({txt.stat().st_size / 1e6:.1f} MB)")
    else:
        print(f"  have {rel(txt)} (use --refresh-mac to update)")

    with txt.open(encoding="utf-8", errors="replace") as f:
        first = f.readline().strip().strip('"')  # NMDP quotes this line
    last_updated = first if first.upper().startswith("LAST") else None
    code_line = re.compile(r"^\*?\s*[A-Z]+\s+\S*\d")
    with txt.open(encoding="utf-8", errors="replace") as f:
        n_codes = sum(1 for line in f if code_line.match(line))
    return {
        "url": MAC_URL,
        "file": rel(txt),
        "last_updated": last_updated,
        "codes": n_codes,
        "sha256": sha256(txt),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--release", default=PINNED_RELEASE, help=f"IMGT/HLA release (default {PINNED_RELEASE})")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--force", action="store_true", help="re-download IMGT/HLA files")
    parser.add_argument("--refresh-mac", action="store_true", help="re-download the MAC table")
    parser.add_argument("--skip-mac", action="store_true", help="IMGT/HLA files only")
    args = parser.parse_args()

    data = args.data_dir.resolve()
    print(f"IMGT/HLA {args.release}")
    manifest = {"imgthla": fetch_imgt(args.release, data, args.force)}

    derived_dir = data / "derived" / args.release
    print(f"Lookup tables -> {rel(derived_dir)}")
    counts = derive_all(data / "imgthla" / args.release, derived_dir)
    manifest["derived"] = {"dir": rel(derived_dir), "rows": counts}
    if not args.skip_mac:
        print("NMDP multiple-allele codes")
        manifest["mac"] = fetch_mac(data, args.refresh_mac)
    manifest["built"] = datetime.now(UTC).isoformat(timespec="seconds")

    out = data / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {rel(out)}")


if __name__ == "__main__":
    main()
