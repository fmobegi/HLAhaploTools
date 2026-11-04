"""Rebuild typed parents from their children and score them against the parents' own typing.

For families in which both parents and at least one child are typed and the pedigree checks find
no problem, the parents are hidden, both together and one at a time, and rebuilt from the
remaining members. Each rebuilt allele is compared with that parent's typing as the tool reads it
(blank second alleles and DRB3/4/5 filled as for any input), at the target resolution:

    correct     the allele is in the parent's typing, in the pairing that matches best
    incorrect   it is not
    unassigned  the children did not determine it (LOCUS*?)
    in place of a blank   the parent's second allele was left blank (read as homozygous) and the
                children show another allele there: a typing gap, not a wrong rebuild

With both parents hidden, which rebuilt parent is the father is arbitrary, so the better of the two
assignments is scored. Families in which the full phasing reports a recombination or a locus that
no inheritance pattern fits are scored apart ("flagged"): without a recombination model, a crossover
near one end of the typed region can make a rebuilt parent wrong there.

This is the larger-cohort counterpart of validate_pedigrees.py, which also tests trios and
parentage on small, fully checked pedigrees.

    uv run python scripts/validate_rebuilt_parents.py FAMILY_TABLE [--resolution 2f] [--loci A,C,B,...] \\
        [--out bench/rebuilt_parents]

Writes OUT/summary.tsv (per scenario, family group and locus group) and OUT/alleles.tsv (every
scored allele, for review), and prints the summary.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import hlahaplotools as h

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_pedigrees import clean  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SERIOUS = {"parent_excluded", "trio_inconsistent", "single_locus_mismatch", "inheritance_mismatch",
           "not_full_sibling", "possible_duplicate"}  # fmt: skip
CLASSICAL = {"A", "B", "C", "DRB1", "DQB1"}


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    text = path.read_text(encoding="utf-8-sig")
    dialect = csv.Sniffer().sniff(text.splitlines()[0], delimiters=",\t;")
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    return list(reader.fieldnames or []), list(reader)


def role(r: dict) -> str:
    m = r["Family_Member"].strip().upper()
    if m in ("F", "FATHER"):
        return "F"
    if m in ("M", "MOTHER"):
        return "M"
    return "C" if re.fullmatch(r"C\d*|CHILD\d*", m) else "?"


def sample_id(r: dict) -> str:
    return (r.get("Sample_ID") or "").strip() or f"{r['FAMILY_ID'].strip()}/{r['Family_Member'].strip()}"


def phase(tmp: Path, header, rows, nom, res, mac, loci):
    p = tmp / "table.csv"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    data = h.read_family_table(p, mac_codes=mac, loci=loci)
    return data, h.phase_families(data, nom, res)


def typed_sides(genotype: str, res: str, nom) -> list[set[str]]:  # res: one locus's resolution
    """'A*01:01+A*02:01' (alternatives '|', alleles '/') -> the two sides as sets of names at `res`."""
    sides: list[set[str]] = [set(), set()]
    for alt in genotype.split("|"):
        for k, side in enumerate(alt.split("+")[:2]):
            for part in side.split("/"):
                if not part:
                    continue
                try:
                    sides[k].update(c for c, _ in nom.translate(part, res))
                except ValueError:
                    sides[k].add(part)  # absent gene (NNNN) or a locus outside the tables
    return sides


def score(rebuilt: list[str], sides: list[set[str]]) -> tuple[int, int, int]:
    """Best pairing of two rebuilt alleles with the two typed sides: (correct, incorrect, unassigned)."""
    best = None
    for a, b in ((0, 1), (1, 0)):
        res = []
        for allele, side in ((rebuilt[a], sides[0]), (rebuilt[b], sides[1])):
            res.append(None if allele.endswith("*?") else allele in side)
        key = (sum(r is True for r in res), -sum(r is False for r in res))
        if best is None or key > best[0]:
            best = (key, res)
    r = best[1]
    return sum(x is True for x in r), sum(x is False for x in r), sum(x is None for x in r)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", type=Path, help="family table (.csv or .tsv)")
    ap.add_argument("--resolution", default="2f")
    ap.add_argument("--loci", help="comma-separated loci (default: every locus in the table)")
    ap.add_argument("--out", type=Path, default=ROOT / "bench" / "rebuilt_parents")
    args = ap.parse_args()

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
    mac = str(ROOT / "data" / "mac" / "alpha.v3.txt")
    header, rows = read_rows(args.table)
    # loci outside the IMGT/HLA tables (MICA, MICB, H, J...) are used as written, cut to two fields
    in_table, uncovered, _ = clean(header, rows, nom, args.resolution)
    loci = args.loci.split(",") if args.loci else in_table
    res = {l: ("asItIs" if l in uncovered else args.resolution) for l in loci}
    by_family = defaultdict(list)
    for r in rows:
        if r.get("FAMILY_ID", "").strip():
            by_family[r["FAMILY_ID"].strip()].append(r)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        data, ref = phase(tmp, header, rows, nom, res, mac, loci)
        ds = data.dataset
        as_read = {sid: ds.sample(i) for i, sid in enumerate(ds.ids)}
        mendel = {m["family"] for m in h.check_mendelian(data, nom, res) if m["locus"]}
        flags = h.audit_families(data, nom, res)
        serious = {f["family"] for f in flags if f["kind"] in SERIOUS}
        flagged = {f["family"] for f in flags if f["kind"] == "recombination"}
        flagged |= {p.family for p in ref if p.inconsistent_loci}

        eligible = []
        for fam, members in by_family.items():
            roles = Counter(role(r) for r in members)
            if roles["F"] == 1 and roles["M"] == 1 and roles["C"] >= 1 and fam not in mendel and fam not in serious:
                eligible.append(fam)
        print(
            f"{len(by_family)} families; {len(eligible)} with both parents and a child typed and no pedigree "
            f"problem ({sum(f in flagged for f in eligible)} of them with recombination or an inconsistent locus)"
        )

        totals: dict[tuple, Counter] = defaultdict(Counter)
        detail = []
        for scenario, hidden in (("both parents hidden", "FM"), ("father hidden", "F"), ("mother hidden", "M")):
            keep = [r for f in eligible for r in by_family[f] if role(r) not in hidden]
            _, phases = phase(tmp, header, keep, nom, res, mac, loci)
            for p in phases:
                parents = {role(r): r for r in by_family[p.family] if role(r) in "FM"}
                halves = {"F": p.haplotypes[0:2], "M": p.haplotypes[2:4]}
                options = [halves] if len(hidden) == 1 else [halves, {"F": halves["M"], "M": halves["F"]}]
                best = None
                for opt in options:
                    rows_, tally = [], Counter()
                    for who in hidden:
                        truth = as_read.get(sample_id(parents[who]), {})
                        # by locus name: the haplotypes hold only the loci the tool phased
                        haps = [{a.split("*")[0]: a for a in hap.split("~")} for hap in opt[who]]
                        for locus in loci:
                            if not truth.get(locus) or locus not in haps[0] or locus not in haps[1]:
                                continue
                            pair = [haps[0][locus], haps[1][locus]]
                            sides = typed_sides(truth[locus], res[locus], nom)
                            c, i, u = score(pair, sides)
                            # a second allele left blank in the file and read as homozygous: an allele the
                            # children show in its place is the blank, not a wrong rebuild
                            raw = parents[who]
                            one_typed = sum(bool((raw.get(f"{locus}_{j}") or "").strip(" *NA")) for j in (1, 2)) == 1
                            blank = i if one_typed and locus not in ("DRB3", "DRB4", "DRB5") else 0
                            i -= blank
                            tally.update({"correct": c, "incorrect": i})
                            rows_.append(
                                [scenario, p.family, who, sample_id(parents[who]), locus, truth[locus],
                                 "+".join(pair), c, i, u, blank]
                            )  # fmt: skip
                    if best is None or (tally["correct"], -tally["incorrect"]) > (
                            best[1]["correct"], -best[1]["incorrect"]):  # fmt: skip
                        best = (rows_, tally)
                for row in best[0]:
                    group = "flagged" if row[1] in flagged else "clean"
                    loci_group = "classical" if row[4] in CLASSICAL else "other"
                    for key in ((scenario, group, "all loci"), (scenario, group, loci_group)):
                        totals[key].update(
                            {"correct": row[7], "incorrect": row[8], "unassigned": row[9], "blank": row[10]}
                        )
                    detail.append(row + [group])

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "alleles.tsv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "scenario",
                "family",
                "parent",
                "sample",
                "locus",
                "typing_as_read",
                "rebuilt",
                "correct",
                "incorrect",
                "unassigned",
                "in_place_of_blank",
                "family_group",
            ]
        )
        w.writerows(detail)
    with (args.out / "summary.tsv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["scenario", "families", "loci", "correct", "incorrect", "unassigned", "in_place_of_blank",
                    "correct_of_assigned"])  # fmt: skip
        for (scenario, group, lg), t in sorted(totals.items()):
            assigned = t["correct"] + t["incorrect"]
            pct = f"{100 * t['correct'] / assigned:.2f}" if assigned else ""
            w.writerow([scenario, group, lg, t["correct"], t["incorrect"], t["unassigned"], t["blank"], pct])
            if lg == "all loci":
                print(
                    f"  {scenario:20} {group:8} correct {t['correct']:6}  incorrect {t['incorrect']:4}  "
                    f"unassigned {t['unassigned']:6}  in place of a blank {t['blank']:3}  "
                    f"correct of assigned {pct}%"
                )
    print(f"Wrote {args.out / 'summary.tsv'} and {args.out / 'alleles.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
