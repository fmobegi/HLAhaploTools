"""Does the EM find the haplotypes that the families show?

Families are phased from inheritance (the gold standard). Then the same people go through the EM as if they
were unrelated, and two questions are asked:

1. Ranking: are the haplotypes the EM puts at the top of its list family haplotypes?
2. Pairs: is each person's most likely EM haplotype pair the pair the family shows?

A family haplotype may leave an allele unresolved (`*?`); it then counts as "any allele", and only the loci
known in the family are compared. Families with Mendelian inconsistencies are left out. Two EM inputs are
tried: everyone (children carry their parents' haplotypes, so the EM sees each one more than once) and the
founders alone (each haplotype seen once or twice, the hardest case). Each locus set is estimated on its own
loci; the gold standard is the 18-locus phasing projected onto those loci.

    uv run python scripts/validate_family_em.py examples/families/*.csv [--json bench/family_em.json]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path

import hlahaplotools as h

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_pedigrees import clean, load_rows  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ROUTINE11 = ["A", "B", "C", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"]
CLASSICAL6 = ["A", "B", "C", "DRB1", "DQB1", "DPB1"]


def write(path: Path, header: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def split(hap: str) -> dict[str, str | None]:
    """Locus -> allele, None where the family left it unresolved."""
    out = {}
    for a in hap.split("~"):
        locus = a.split("*", 1)[0]
        out[locus] = None if a.endswith("*?") else a
    return out


def project(hap: dict[str, str | None], loci: list[str]) -> tuple:
    return tuple(hap.get(l) for l in loci)


def compatible(em: tuple, gold: tuple) -> bool:
    return all(g is None or e == g for e, g in zip(em, gold, strict=True))


def top_share(in_gold: list[bool], k: int) -> float:
    """Share of the k most frequent EM haplotypes that are family haplotypes."""
    return sum(in_gold[:k]) / min(k, len(in_gold)) if in_gold else float("nan")


def pair_score(em: tuple[tuple, tuple], gold: tuple[tuple, tuple]) -> tuple[int, int]:
    """(loci where the EM pair agrees with the family pair, loci the family knows on both haplotypes),
    in the better of the two orientations."""
    best = (0, 0)
    for e1, e2 in (em, em[::-1]):
        agree = known = 0
        for k in range(len(gold[0])):
            g1, g2 = gold[0][k], gold[1][k]
            if g1 is None or g2 is None:
                continue
            known += 1
            agree += e1[k] == g1 and e2[k] == g2
        best = max(best, (agree, known))
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tables", nargs="+", type=Path)
    ap.add_argument("--resolution", default="2f")
    ap.add_argument("--starts", type=int, default=10)
    ap.add_argument("--json", type=Path, help="also write all results as JSON here")
    args = ap.parse_args()

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
    mac = str(ROOT / "data" / "mac" / "alpha.v3.txt")
    header, rows = load_rows(args.tables)
    loci, uncovered, _ = clean(header, rows, nom, args.resolution)
    res_all = {l: ("asItIs" if l in uncovered else args.resolution) for l in loci}

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        write(tmp / "all.csv", header, rows)
        data = h.read_family_table(tmp / "all.csv", mac_codes=mac)
        phases = h.phase_families(data, nom, res_all)
        order = list(data.dataset.loci)

        # ---- gold standard: founders' haplotypes and each person's pair, from clean families
        bad = {p.family for p in phases if p.inconsistent_loci or p.status == "failed"}
        gold_haps: list[dict] = []
        gold_pairs: dict[str, tuple[dict, dict]] = {}
        for p in phases:
            if p.family in bad:
                continue
            haps = [split(x) for x in p.haplotypes]
            gold_haps += haps
            labels = dict(zip("ABCD", haps, strict=True))
            for k in range(2):
                pid = p.parent_id(k)
                if any(r for r in rows if f"{r['FAMILY_ID']}/{r['Family_Member']}" == pid):
                    gold_pairs[pid] = (haps[2 * k], haps[2 * k + 1])
            for cid, pat, mat in p.children:
                if pat and mat:
                    gold_pairs[cid] = (labels[pat], labels[mat])
        keep = [r for r in rows if r["FAMILY_ID"] not in bad]
        founders = [r for r in keep if r["Family_Member"].strip().upper() in ("F", "M", "FATHER", "MOTHER")]
        print(f"Families: {len(phases)}, left out (Mendelian inconsistencies): {', '.join(sorted(bad)) or '-'}")
        print(f"People: {len(keep)} kept, of whom {len(founders)} founders; family pairs known for {len(gold_pairs)}")

        sets = {
            "18 loci (all typed)": list(order),
            "11 routinely typed loci": [l for l in order if l in ROUTINE11],
            "6 classical loci": [l for l in order if l in CLASSICAL6],
        }
        report = {"families_left_out": sorted(bad), "people": len(keep), "founders": len(founders), "results": []}
        for set_name, set_loci in sets.items():
            gold = {project(g, set_loci) for g in gold_haps}
            gold_full = {g for g in gold if None not in g}
            for input_name, people in (("everyone", keep), ("founders only", founders)):
                cols = ["FAMILY_ID", "Family_Member"] + [c for c in header if c.rsplit("_", 1)[0] in set_loci]
                path = tmp / "subset.csv"
                write(path, cols, people)
                sub = h.read_family_table(path, mac_codes=mac)
                res = {l: res_all[l] for l in set_loci}
                est = h.estimate(sub.dataset, nom, res, starts=args.starts, epsilon=1e-6)
                sub_order = list(sub.dataset.loci)
                assert sub_order == set_loci, (sub_order, set_loci)
                ranked = [(tuple(split(hp).get(l) for l in set_loci), f) for hp, f in est.frequencies(cut=0.0)]
                ranked.sort(key=lambda x: -x[1])
                in_gold = [any(compatible(hp, g) for g in gold) for hp, _ in ranked]
                n = len(gold)
                mass = sum(f for (hp, f), ok in zip(ranked, in_gold, strict=True) if ok)
                found_full = sum(1 for g in gold_full if any(hp == g for hp, _ in ranked))
                # Pairs: each kept person with a known family pair and an EM pair.
                counts = Counter()
                for sid, h1, h2, _prob in est.best_diplotypes():
                    if sid not in gold_pairs:
                        continue
                    g = tuple(project(x, set_loci) for x in gold_pairs[sid])
                    e = (tuple(split(h1).get(l) for l in set_loci), tuple(split(h2).get(l) for l in set_loci))
                    agree, known = pair_score(e, g)
                    counts["people"] += 1
                    counts["exact"] += agree == known
                    counts["loci_agree"] += agree
                    counts["loci_known"] += known
                row = {
                    "locus_set": set_name, "input": input_name, "people_in_em": est.n_samples,
                    "em_haplotypes": len(ranked), "family_haplotypes": n,
                    "top10_in_family": top_share(in_gold, 10), "topN_in_family": top_share(in_gold, n),
                    "em_frequency_on_family_haplotypes": mass,
                    "complete_family_haplotypes_found": f"{found_full}/{len(gold_full)}",
                    "pairs_compared": counts["people"],
                    "pairs_exact": counts["exact"],
                    "locus_agreement": counts["loci_agree"] / max(counts["loci_known"], 1),
                }  # fmt: skip
                report["results"].append(row)

    print()
    print(
        "| Locus set | EM input | People | Family haplotypes | Top 10 in family | Top N in family | "
        "EM frequency on family haplotypes | Pairs exact | Loci agreeing |"
    )
    print("|---|---|---|---|---|---|---|---|---|")
    for r in report["results"]:
        print(f"| {r['locus_set']} | {r['input']} | {r['people_in_em']} | {r['family_haplotypes']} | "
              f"{r['top10_in_family']:.0%} | {r['topN_in_family']:.0%} | "
              f"{r['em_frequency_on_family_haplotypes']:.1%} | "
              f"{r['pairs_exact']}/{r['pairs_compared']} | {r['locus_agreement']:.1%} |")  # fmt: skip
    print(
        "\nTop N: the N most frequent EM haplotypes, N = the number of distinct family haplotypes. A family "
        "haplotype with an unresolved allele matches any allele there."
    )
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
