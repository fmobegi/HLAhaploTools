"""How well does audit_families find planted pedigree errors?

Simulates families with known problems (scripts/simulate.py --errors), runs
the audit, and reports for each kind of problem how often it is flagged at all
and how often with the expected flag, plus how often clean children are
flagged (false alarms). Results are split by whether both parents are typed.

    uv run python scripts/validate_audit.py --families 2000 --seed 1
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import hlahaplotools as h

ROOT = Path(__file__).resolve().parent.parent

EXPECTED = {
    "other_father": {"parent_excluded", "not_full_sibling", "inheritance_mismatch"},
    "duplicate_of_mother": {"possible_duplicate"},
    "recombinant": {"recombination"},
    "typo": {"single_locus_mismatch", "trio_inconsistent", "inheritance_mismatch", "recombination"},
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", type=int, default=2000)
    ap.add_argument("--children", default="2-4")
    ap.add_argument("--missing-parent", type=float, default=0.15)
    ap.add_argument("--loci", default="A,C,B,DRB1,DQB1,DPB1")
    ap.add_argument("--errors", default="other_father:0.04,duplicate_of_mother:0.02,recombinant:0.04,typo:0.04")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "simulate.py"),
                str(tmp),
                "--format",
                "families",
                "--samples",
                str(args.families),
                "--children",
                args.children,
                "--missing-parent",
                str(args.missing_parent),
                "--loci",
                args.loci,
                "--errors",
                args.errors,
                "--seed",
                str(args.seed),
                "--haplotypes",
                "1500",
                "--alleles",
                "40",
            ],
            check=True,
            capture_output=True,
        )
        planted = {}
        for line in (tmp / "truth_errors.tsv").read_text().splitlines():
            _, member, kind = line.split("\t")
            planted[member] = kind
        data = h.read_family_table(tmp / "families.csv")
        flags = h.audit_families(data, nom, "asItIs")

    fams = data.families()
    both_parents = {f["family"]: f["father"] is not None and f["mother"] is not None for f in fams}
    children = [(f["family"], c) for f in fams for c in f["children"]]
    by_member = defaultdict(set)
    for fl in flags:
        if fl["sample"]:
            by_member[fl["sample"]].add(fl["kind"])
        else:  # family-level flag naming candidate children
            for _, c in children:
                if c.split("/")[0] == fl["family"] and f"{c}," in fl["message"] + ",":
                    by_member[c].add(fl["kind"] + " (one of several)")

    rows = Counter()
    for fam, c in children:
        kind = planted.get(c, "none")
        setting = "both parents" if both_parents[fam] else "a parent missing"
        got = by_member.get(c, set())
        rows[(kind, setting, "n")] += 1
        rows[(kind, setting, "flagged")] += bool(got)
        if kind != "none":
            rows[(kind, setting, "expected")] += bool({g.split(" (")[0] for g in got} & EXPECTED[kind])
        else:
            for g in got:
                rows[("none", setting, f"false:{g}")] += 1

    print(
        f"{args.families} families, {len(children)} children, {len(flags)} flags; "
        f"loci {args.loci}; {args.missing_parent:.0%} of parents untyped\n"
    )
    print(f"{'planted problem':22} {'parents':17} {'children':>8} {'flagged':>9} {'expected flag':>14}")
    report = []
    for kind in [*EXPECTED, "none"]:
        for setting in ("both parents", "a parent missing"):
            n = rows[(kind, setting, "n")]
            if not n:
                continue
            flagged = rows[(kind, setting, "flagged")]
            line = f"{kind:22} {setting:17} {n:8} {flagged / n:9.1%}"
            if kind != "none":
                line += f" {rows[(kind, setting, 'expected')] / n:14.1%}"
            else:
                line += "   (false alarms)"
            print(line)
            report.append(
                {
                    "planted": kind,
                    "parents": setting,
                    "n": n,
                    "flagged": flagged,
                    "expected": rows[(kind, setting, "expected")],
                }
            )
    false = {k[2][6:]: v for k, v in rows.items() if k[0] == "none" and k[2].startswith("false:")}
    if false:
        print("\nflags on clean children, by kind:", dict(sorted(false.items(), key=lambda x: -x[1])))
    if args.json:
        args.json.write_text(json.dumps({"rows": report, "false_alarms": false}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
