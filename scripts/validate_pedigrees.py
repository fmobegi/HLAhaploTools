"""Test family phasing on real pedigrees by leaving members out.

For every family the reference is the phasing with *all* members. Then:

1. Trios: both parents + one child at a time. Are the parental haplotypes the
   same as the reference? (Only loci determined in both are scored.)
2. Parentage: does each true child pass the Mendel check with its parents, and
   does every child from another family fail?
3. One parent missing: drop the father (then the mother), keep 2+ children.
   Is the rebuilt parent the same as the reference, and are its alleles part of
   that parent's actual typing?

Usage (from the project root):

    uv run python scripts/validate_pedigrees.py families1.csv families2.csv ... [--resolution 2f]

Alleles at loci the nomenclature tables do not cover (H, J, MICA, MICB) are cut
to 2 fields and used as written; single-digit first fields (Excel dropping a
leading zero, `1:01:14`) are restored to `01:01:14`. Both are counted.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import hlahaplotools as h

ROOT = Path(__file__).resolve().parent.parent
FIELD_FIX = re.compile(r"^(\*?)(\d)(:)")


def load_rows(paths: list[Path]) -> tuple[list[str], list[dict[str, str]]]:
    """All rows of all files, with the union of columns (files may order them differently)."""
    header: list[str] = []
    rows: list[dict[str, str]] = []
    for p in paths:
        with p.open(newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            for col in reader.fieldnames or []:
                if col not in header:
                    header.append(col)
            rows.extend(reader)
    return header, rows


def clean(header, rows, nom, resolution):
    """Restore dropped leading zeros; cut uncovered loci to 2 fields. Returns fixes per kind."""
    loci = sorted({c.rsplit("_", 1)[0] for c in header if re.search(r"_[12]$", c)})
    uncovered = set()
    for locus in loci:
        sample = next((r[f"{locus}_1"] for r in rows if r.get(f"{locus}_1", "").strip("* ") not in ("", "NA")), None)
        if sample is None:
            continue
        name = sample if "*" in sample and not sample.startswith("*") else f"{locus}*{sample.lstrip('*')}"
        try:
            nom.translate(FIELD_FIX.sub(r"\g<1>0\2\3", name.split("*", 1)[0] + "*" + name.split("*", 1)[1]), resolution)
        except ValueError:
            uncovered.add(locus)
    fixes = Counter()
    for r in rows:
        for col, v in r.items():
            if not re.search(r"_[12]$", col) or v is None:
                continue
            v = v.strip()
            if FIELD_FIX.match(v):
                v = FIELD_FIX.sub(r"\g<1>0\2\3", v)
                fixes["leading zero restored"] += 1
            locus = col.rsplit("_", 1)[0]
            if locus in uncovered and v not in ("", "NA"):
                prefix, _, typing = v.rpartition("*")
                cut = ":".join(typing.split(":")[:2])
                if cut != typing:
                    fixes[f"{locus} cut to 2 fields"] += 1
                v = f"{prefix}*{cut}" if prefix or v.startswith("*") else cut
            r[col] = v
    return loci, uncovered, fixes


def write(path: Path, header, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def phase_rows(tmp: Path, header, rows, nom, res, mac):
    p = tmp / "subset.csv"
    write(p, header, rows)
    data = h.read_family_table(p, mac_codes=mac)
    return h.phase_families(data, nom, res), h.check_mendelian(data, nom, res), data


def alleles(hap: str) -> list[str]:
    return hap.split("~")


def score_parent(test: list[str], ref: list[str]) -> dict:
    """Compare a parent's two haplotypes (test) with the reference, in the better pairing.

    Per locus and haplotype: correct, wrong, unknown in test, or not scorable
    (unknown in the reference).
    """
    best = None
    for t1, t2 in (test, test[::-1]):
        c = Counter()
        for th, rh in ((t1, ref[0]), (t2, ref[1])):
            for ta, ra in zip(alleles(th), alleles(rh), strict=True):
                if ra.endswith("*?"):
                    c["not scorable"] += 1
                elif ta.endswith("*?"):
                    c["unknown"] += 1
                elif ta == ra:
                    c["correct"] += 1
                else:
                    c["wrong"] += 1
        if best is None or (c["correct"], -c["wrong"]) > (best["correct"], -best["wrong"]):
            best = c
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tables", nargs="+", type=Path)
    ap.add_argument("--resolution", default="2f")
    ap.add_argument("--json", type=Path, help="also write all results as JSON here")
    args = ap.parse_args()

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
    mac = str(ROOT / "data" / "mac" / "alpha.v3.txt")
    header, rows = load_rows(args.tables)
    loci, uncovered, fixes = clean(header, rows, nom, args.resolution)
    res = {l: ("asItIs" if l in uncovered else args.resolution) for l in loci}
    by_family = defaultdict(list)
    for r in rows:
        by_family[r["FAMILY_ID"]].append(r)
    role = lambda r: r["Family_Member"].strip().upper()  # noqa: E731
    is_parent = lambda r, p: role(r) in ({"F", "FATHER"} if p == "F" else {"M", "MOTHER"})  # noqa: E731
    is_child = lambda r: re.fullmatch(r"C\d*|CHILD\d*", role(r)) is not None  # noqa: E731

    print(f"Loci ({len(loci)}): {', '.join(loci)}")
    print(f"Resolution {args.resolution}; used as written (not in IMGT tables): {', '.join(sorted(uncovered)) or '-'}")
    for k, n in fixes.items():
        print(f"  fixed: {k} x{n}")
    report = {"loci": loci, "uncovered": sorted(uncovered), "fixes": dict(fixes)}

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # ---- reference: all members
        ref_phases, ref_mendel, data = phase_rows(tmp, header, rows, nom, res, mac)
        ref = {p.family: p for p in ref_phases}
        print(f"\nREFERENCE (all members): {len(ref_phases)} families, {len(data)} people")
        for p in ref_phases:
            n_unknown = sum(a.endswith("*?") for hap in p.haplotypes for a in alleles(hap)) if p.haplotypes[0] else None
            print(
                f"  {p.family:10} {p.status:8} children={len(p.children)} "
                f"inconsistent={','.join(p.inconsistent_loci) or '-'} unknown alleles={n_unknown}"
            )
        print(f"  Mendel problems: {len(ref_mendel)}")
        for m in ref_mendel:
            print(f"    {m['family']} {m['sample'] or ''} {m['locus'] or ''}: {m['message']}")
        report["reference"] = [
            {
                "family": p.family,
                "status": p.status,
                "haplotypes": p.haplotypes,
                "inconsistent": p.inconsistent_loci,
                "notes": p.notes,
            }
            for p in ref_phases
        ]
        report["reference_mendel"] = ref_mendel

        # ---- 1. trios
        trio_totals, trio_rows = Counter(), []
        for fam, members in by_family.items():
            parents = [r for r in members if is_parent(r, "F") or is_parent(r, "M")]
            if len(parents) < 2 or ref[fam].status == "failed":
                continue
            for child in (r for r in members if is_child(r)):
                (p,), _, _ = phase_rows(tmp, header, parents + [child], nom, res, mac)
                for i, who in enumerate("FM"):
                    s = score_parent(p.haplotypes[2 * i : 2 * i + 2], ref[fam].haplotypes[2 * i : 2 * i + 2])
                    trio_totals.update(s)
                    trio_rows.append({"family": fam, "child": child["Family_Member"], "parent": who, **s})
        print("\n1. TRIOS (both parents + one child) vs. reference, per haplotype allele")
        _print_scores(trio_totals)
        worst = sorted(trio_rows, key=lambda r: -r.get("wrong", 0))[:5]
        for r in worst:
            if r.get("wrong"):
                print(f"   wrong: {r['family']} child {r['child']} parent {r['parent']}: {r['wrong']} alleles")
        report["trios"] = {"totals": dict(trio_totals), "rows": trio_rows}

        # ---- 2. parentage
        print("\n2. PARENTAGE: Mendel check of every child against every complete parent pair")
        pairs = {f: [r for r in m if is_parent(r, "F") or is_parent(r, "M")] for f, m in by_family.items()}
        pairs = {f: ps for f, ps in pairs.items() if len(ps) == 2}
        children = [(f, r) for f, m in by_family.items() for r in m if is_child(r)]
        outcome = Counter()
        mismatch_counts = defaultdict(list)
        for fam, parents in pairs.items():
            for cfam, child in children:
                c = dict(child, FAMILY_ID=fam)
                _, mendel, _ = phase_rows(tmp, header, parents + [c], nom, res, mac)
                n = sum(1 for m in mendel if m["locus"])
                true_child = cfam == fam
                mismatch_counts["true" if true_child else "false"].append(n)
                outcome[("true child" if true_child else "other family", "pass" if n == 0 else "fail")] += 1
        for (kind, verdict), n in sorted(outcome.items()):
            print(f"   {kind:12} {verdict}: {n}")
        t, f = mismatch_counts["true"], mismatch_counts["false"]
        print(
            f"   inconsistent loci per child: true children {sorted(t)}; other families min {min(f)} "
            f"median {sorted(f)[len(f) // 2]} max {max(f)}"
        )
        report["parentage"] = {f"{k[0]} {k[1]}": n for k, n in outcome.items()} | {
            "true_child_mismatches": t,
            "other_family_mismatches": f,
        }

        # ---- 3. one parent missing, 2+ children
        print("\n3. ONE PARENT MISSING (2+ children) vs. reference")
        miss_totals = {"rebuilt": Counter(), "known": Counter()}
        geno = Counter()
        miss_rows = []
        for fam, members in by_family.items():
            kids = [r for r in members if is_child(r)]
            if fam not in pairs or len(kids) < 2:
                continue
            for drop in "FM":
                keep = [r for r in members if not is_parent(r, drop)]
                for n_kids in sorted({2, len(kids)}):
                    for subset in itertools.combinations(kids, n_kids):
                        rows_used = [r for r in keep if not is_child(r)] + list(subset)
                        (p,), _, _ = phase_rows(tmp, header, rows_used, nom, res, mac)
                        for i, who in enumerate("FM"):
                            kind = "rebuilt" if who == drop else "known"
                            s = score_parent(p.haplotypes[2 * i : 2 * i + 2], ref[fam].haplotypes[2 * i : 2 * i + 2])
                            miss_totals[kind].update(s)
                            miss_rows.append(
                                {"family": fam, "dropped": drop, "children": n_kids, "parent": who, "kind": kind, **s}
                            )
                            if kind == "rebuilt":
                                actual = next(r for r in members if is_parent(r, drop))
                                geno.update(_genotype_check(p.haplotypes[2 * i : 2 * i + 2], actual, res, nom))
        for kind in ("rebuilt", "known"):
            print(f"   {kind} parent:")
            _print_scores(miss_totals[kind], indent="     ")
        print("   rebuilt parent's alleles vs. that parent's own typing (same resolution):")
        total = sum(geno.values())
        for k in ("in typing", "not in typing", "parent untyped here"):
            print(f"     {k:20} {geno[k]:5} ({geno[k] / total:.1%})")
        report["one_parent_missing"] = {
            "rebuilt": dict(miss_totals["rebuilt"]),
            "known": dict(miss_totals["known"]),
            "rebuilt_vs_own_typing": dict(geno),
            "rows": miss_rows,
        }

        # ---- children-only families
        orphan = [p for p in ref_phases if p.father is None and p.mother is None]
        if orphan:
            print("\n4. FAMILIES WITHOUT PARENTS")
            for p in orphan:
                n = sum(not a.endswith("*?") for hap in p.haplotypes for a in alleles(hap))
                tot = sum(1 for hap in p.haplotypes for _ in alleles(hap))
                print(f"   {p.family}: {p.status}, {n}/{tot} parental alleles recovered; {'; '.join(p.notes)}")

    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + "\n")
    return 0


def _genotype_check(haps: list[str], actual: dict, res: dict, nom) -> Counter:
    """Is each recovered allele of a rebuilt parent part of that parent's real typing?"""
    out = Counter()
    for a1, a2 in zip(alleles(haps[0]), alleles(haps[1]), strict=True):
        locus = a1.split("*")[0]
        typed = set()
        for col in (f"{locus}_1", f"{locus}_2"):
            v = (actual.get(col) or "").strip()
            if v in ("", "NA"):
                continue
            name = v if "*" in v and not v.startswith("*") else f"{locus}*{v.lstrip('*')}"
            for part in name.split("/"):
                part = part if part.startswith(locus + "*") else f"{locus}*{part.split('*')[-1]}"
                try:
                    typed.update(c for c, _ in nom.translate(part, res[locus]))
                except ValueError:
                    typed.add(part)
        for a in (a1, a2):
            if a.endswith("*?"):
                continue
            out["parent untyped here" if not typed else ("in typing" if a in typed else "not in typing")] += 1
    return out


def _print_scores(c: Counter, indent: str = "   ") -> None:
    scored = c["correct"] + c["wrong"]
    called = scored + c["unknown"]
    print(
        f"{indent}correct {c['correct']}, wrong {c['wrong']}, left unknown {c['unknown']}, "
        f"not scorable {c['not scorable']}"
    )
    if scored:
        print(
            f"{indent}accuracy of called alleles: {c['correct'] / scored:.2%}; "
            f"called {scored / called:.1%} of scorable alleles"
        )


if __name__ == "__main__":
    raise SystemExit(main())
