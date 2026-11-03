"""Update allele names that IMGT/HLA has since deleted or renamed.

IMGT/HLA records every retired name in ``Deleted_alleles.txt`` with a note such
as "Suffix changed becoming C*04:09L" or "Sequence shown to be identical to
B*47:01:01:03". :func:`load_renames` turns those notes into an old -> current
name table (keeping only targets present in the release's ``Allelelist.txt``),
and :func:`rename_file` rewrites a typing file with the current names, so the
file can then be estimated at any resolution, or as typed.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

# The name a retired allele now goes by, from the wording IMGT/HLA uses.
_TARGET = re.compile(
    r"(?:becoming|identical to(?: the)?(?: corrected)?(?: sequence of)?|renamed(?: to)?)\s+"
    r"([A-Z][A-Z0-9]*\*\d+(?::\d+)+[A-Z]?)"
)
# Locus column names in family tables: A_1, A.1, A1, A_ALLELE_1, HLA-A_1.
_COLUMN = re.compile(r"^(?:HLA-)?([A-Za-z0-9]+?)[_.]?(?:ALLELE_?)?[12]$", re.IGNORECASE)


def _rows(path: Path):
    """CSV rows of an IMGT/HLA list file, skipping its '#' header lines."""
    with path.open(encoding="utf-8") as f:
        yield from csv.reader(line for line in f if not line.startswith("#"))


def load_renames(imgt_dir: Path) -> dict[str, str]:
    """Old name -> current name, from ``Deleted_alleles.txt`` in ``imgt_dir``."""
    current = {row[1] for row in _rows(imgt_dir / "Allelelist.txt") if len(row) > 1}
    direct: dict[str, str] = {}
    for row in _rows(imgt_dir / "Deleted_alleles.txt"):
        if len(row) < 3 or row[1] == "Allele":
            continue
        m = _TARGET.search(",".join(row[2:]))
        if m:
            direct[row[1]] = m.group(1)
    out = {}
    for old, new in direct.items():
        for _ in range(10):  # follow chains: a target may itself have been renamed later
            if new in current or new not in direct:
                break
            new = direct[new]
        if new in current and new != old:
            out[old] = new
    return out


def _update(cell: str, locus: str | None, renames: dict[str, str], seen: Counter) -> str:
    """Rename every allele in one cell, keeping the cell's own style (prefix or not)."""

    def one(part: str) -> str:
        if not part:
            return part
        if "*" in part and not part.startswith("*"):
            name, bare = part, None
        elif locus:
            bare = part.lstrip("*")
            name = f"{locus}*{bare}"
        else:
            return part
        new = renames.get(name)
        if new is None:
            return part
        seen[(name, new)] += 1
        if bare is None:
            return new
        return ("*" if part.startswith("*") else "") + new.split("*", 1)[1]

    # Split on GL-string operators and keep them.
    return "".join(one(t) if t not in "|+/~^" else t for t in re.split(r"([|+/~^])", cell))


def rename_file(src: Path, dest: Path, fmt: str, renames: dict[str, str]) -> Counter:
    """Write ``src`` to ``dest`` with current allele names. Returns (old, new) -> count.

    ``fmt`` is ``mac`` (tab-separated, header naming each locus twice), ``glsc``
    (ID then GL strings) or ``families`` (CSV/TSV with LOCUS_1/LOCUS_2 columns).
    Everything else in the file (IDs, separators, empty cells) is kept as it is.
    """
    seen: Counter = Counter()
    text = src.read_text(encoding="utf-8-sig", newline="")  # keep \r\n, so the file's line endings are kept
    newline = "\r\n" if "\r\n" in text else "\n"
    lines = text.splitlines()
    out: list[str] = []
    if fmt == "glsc":
        for line in lines:
            parts = re.split(r"(\s+)", line)
            # parts[0] is the ID; GL strings follow, separated by whitespace.
            out.append(parts[0] + "".join(p if p.isspace() else _update(p, None, renames, seen) for p in parts[1:]))
    else:
        if not lines:
            dest.write_text("", encoding="utf-8")
            return seen
        sep = "\t" if fmt == "mac" or "\t" in lines[0] else (";" if lines[0].count(";") > lines[0].count(",") else ",")
        header = lines[0].split(sep)
        if fmt == "mac":
            loci = [None] + [h.strip() or None for h in header[1:]]
        else:
            loci = []
            for h in header:
                m = _COLUMN.match(h.strip().strip('"'))
                loci.append(m.group(1).upper() if m else None)
        out.append(lines[0])
        for line in lines[1:]:
            cells = line.split(sep)
            for i, cell in enumerate(cells):
                locus = loci[i] if i < len(loci) else None
                if locus and cell.strip():
                    quoted = cell.startswith('"') and cell.endswith('"') and len(cell) > 1
                    inner = cell[1:-1] if quoted else cell
                    new = _update(inner, locus, renames, seen)
                    cells[i] = f'"{new}"' if quoted else new
            out.append(sep.join(cells))
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(newline.join(out) + (newline if text.endswith(("\n", "\r")) else ""), encoding="utf-8", newline="")
    return seen
