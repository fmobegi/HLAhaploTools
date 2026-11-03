"""Hapl-o-Mat's genotype list (WRITE_GENOTYPES output, READ input).

Each line is one resolved multilocus genotype of a person:

    id<TAB>type<TAB>weight<TAB>A*01:01g+A*02:01g^B*08:01g+B*07:02g

with the weights of a person's lines summing to 1. `write_genotypes` writes it from a
dataset; `read_genotypes` reads it back into a dataset, to be estimated with resolution
`asItIs` (the codes are already resolved), as Hapl-o-Mat's READ input does.
"""

from __future__ import annotations

import math
import tempfile
from fractions import Fraction
from pathlib import Path

from hlahaplotools._core import read_glsc, resolved_genotypes

MAX_COPIES = 1000  # most genotypes one locus may list once weights are written as repeats


def write_genotypes(path, dataset, nomenclature, resolutions, min_genotype_frequency=1e-5):
    """Write the genotype list; returns the issues (people left out)."""
    rows, issues = resolved_genotypes(dataset, nomenclature, resolutions, min_genotype_frequency=min_genotype_frequency)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for sid, reports in rows:
            for kind, code, w in reports:
                f.write(f"{sid}\t{kind}\t{w:.14g}\t{code}\n")
    return issues


def _locus_cell(genotype: str) -> str | None:
    """'A*01:01+A*?' -> 'A*01:01' (one allele typed); 'A*?+A*?' -> None (untyped)."""
    known = [a for a in genotype.split("+") if not a.endswith("*?")]
    if not known:
        return None
    return "+".join(known) if len(known) == 2 else known[0]


def read_genotypes(path, loci=None):
    """Read a genotype list into a Dataset (and the lines that could not be used).

    The dataset model keeps alternatives per locus, so a person's lines must be all
    combinations of per-locus genotypes with weights that multiply, as Hapl-o-Mat writes
    them; weights are kept exactly as repeats of an alternative. Returns
    (dataset, issues) with issues as (line, id, reason).
    """
    people: dict[str, list[tuple[int, float, list[str]]]] = {}
    issues = []
    order = None
    for n, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        parts = line.split()
        if not parts:
            continue
        if len(parts) != 4:
            issues.append((n, parts[0], "expected id, type, weight and genotype"))
            continue
        sid, _, w, code = parts
        genos = code.split("^")
        names = [g.split("*", 1)[0] for g in genos]
        if order is None:
            order = names
        if names != order:
            issues.append((n, sid, f"loci {names} differ from the first line's {order}"))
            continue
        people.setdefault(sid, []).append((n, float(w), genos))
    loci = list(loci or order or [])
    keep = [order.index(l) for l in loci if order and l in order]
    lines = []
    for sid, rows in people.items():
        marg: list[dict[str, float]] = [dict() for _ in keep]
        for _, w, genos in rows:
            for m, k in zip(marg, keep, strict=True):
                m[genos[k]] = m.get(genos[k], 0.0) + w
        total = sum(w for _, w, _ in rows)
        expected = math.prod(len(m) for m in marg)
        combos = {tuple(g[k] for k in keep) for _, _, g in rows}
        prod_ok = len(combos) == expected and all(
            math.isclose(w / total, math.prod(m[g[k]] / total for m, k in zip(marg, keep, strict=True)), rel_tol=1e-6)
            for _, w, g in rows
        ) if len(rows) == len(combos) else False  # fmt: skip
        if not prod_ok:
            issues.append((rows[0][0], sid, "genotype list is not a product of per-locus lists; not supported"))
            continue
        cells = []
        bad = None
        for locus, m in zip(loci, marg, strict=True):
            fr = {g: Fraction(w / total).limit_denominator(MAX_COPIES) for g, w in m.items()}
            if any(not math.isclose(float(f), m[g] / total, rel_tol=1e-6) for g, f in fr.items()):
                bad = f"weights at {locus} cannot be written as repeats"
                break
            den = math.lcm(*(f.denominator for f in fr.values()))
            if den > MAX_COPIES:
                bad = f"weights at {locus} need more than {MAX_COPIES} repeats"
                break
            alts = [c for g, f in fr.items() for c in [_locus_cell(g)] * int(f * den)]
            if any(a is None for a in alts):
                if all(a is None for a in alts):
                    continue  # untyped locus
                bad = f"{locus} is untyped in some genotypes only"
                break
            cells.append("|".join(alts))
        if bad:
            issues.append((rows[0][0], sid, bad))
            continue
        lines.append(f"{sid}\t" + "\t".join(cells))
    with tempfile.TemporaryDirectory() as tmp:
        glc = Path(tmp) / "genotypes.glc"
        glc.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        # Already resolved: a blank copy stays unknown, and DRB3/4/5 are as written.
        ds = read_glsc(str(glc), loci=loci, half_typed="unknown", drb345=False)
    return ds, issues + [(l, i, r) for l, i, r in ds.issues]


__all__ = ["read_genotypes", "write_genotypes"]
