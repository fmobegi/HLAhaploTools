"""Test the EM of every benchmarked tool against haplotypes known from family segregation.

Families give each person's two haplotypes from inheritance, without any EM. Founders (unrelated
people) are taken from the families, every tool estimates haplotype frequencies from them with
`benchmark.py run`, and each founder's possible haplotype pairs are then ranked by each tool's
frequencies: does the pedigree-derived pair come first?

Three steps:

    # 1. one set of unrelated founders per family: the typed parents, or the first child if neither
    #    parent is typed (optionally leaving out families with pedigree flags)
    uv run python scripts/validate_em_segregation.py founders FAMILY_TABLE -o bench/seg/founders.dat \\
        [--exclude-flags results/families/flags.tsv]

    # 2. the same input for every tool (and haplo.stats on the same cases)
    uv run python scripts/benchmark.py run --out bench/seg --real bench/seg/founders.dat --sizes "" \\
        --real-sets B~C,A~C~B,A~C~B~DRB1~DQB1
    uv run python scripts/bench_haplostats.py --out bench/seg/haplostats --cases bench/seg/cases/real_...

    # 3. score: the family phasing of the same table gives the true pairs
    hlahaplotools families FAMILY_TABLE --loci ... -o results/families
    uv run python scripts/validate_em_segregation.py score --bench bench/seg --families results/families

Only founders heterozygous at two or more loci of a set are scored (otherwise there is one possible
pair). Loci the family run left unphased (LOCUS*?) or marked ambiguous or inconsistent are not used.
Pairs tied at the top share the first rank. Frequencies are in-sample: each founder is part of the
data the frequencies were estimated from. Writes OUT/segregation_check.tsv (per locus set and tool)
and OUT/segregation_people.tsv (per founder and tool).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import re
from collections import defaultdict
from pathlib import Path

FLAG_KINDS = {"parent_excluded", "possible_duplicate"}


def read_rows(path: Path) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8-sig")
    dialect = csv.Sniffer().sniff(text.splitlines()[0], delimiters=",\t;")
    return list(csv.DictReader(text.splitlines(), dialect=dialect))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ------------------------------------------------------------------ founders


def cmd_founders(args) -> int:
    rows = read_rows(args.table)
    loci = list(dict.fromkeys(c.rsplit("_", 1)[0] for c in rows[0] if re.search(r"_[12]$", c)))
    flagged = set()
    for p in args.exclude_flags or []:
        flagged |= {r["family"] for r in read_tsv(p) if r["kind"] in FLAG_KINDS}
    fams = defaultdict(list)
    for r in rows:
        if r.get("FAMILY_ID", "").strip():
            fams[r["FAMILY_ID"].strip()].append(r)

    def role(r):
        m = r["Family_Member"].strip().upper()
        return "P" if m in ("F", "M", "FATHER", "MOTHER") else "C" if re.fullmatch(r"C\d*|CHILD\d*", m) else "?"

    def sid(r):
        return (r.get("Sample_ID") or "").strip() or f"{r['FAMILY_ID'].strip()}/{r['Family_Member'].strip()}"

    used, seen, skipped = [], set(), 0
    for fam, members in fams.items():
        if fam in flagged:
            skipped += 1
            continue
        picks = [m for m in members if role(m) == "P"] or [m for m in members if role(m) == "C"][:1]
        for m in picks:
            if sid(m) not in seen:
                seen.add(sid(m))
                used.append(m)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["id"] + [l for l in loci for _ in (0, 1)])
        for m in used:
            w.writerow([sid(m)] + [(m.get(f"{l}_{k}") or "").split("*")[-1].strip() for l in loci for k in (1, 2)])
    kept = len(fams) - skipped
    print(f"{len(used)} founders from {kept} families ({skipped} flagged families left out) -> {args.output}")
    return 0


# ------------------------------------------------------------------ score


def canon(alleles) -> str:
    return "~".join(sorted(alleles, key=lambda a: a.split("*")[0]))


def bare(hap: str) -> str:
    """A haplotype with expression suffixes (N, L, S, Q...) removed from every allele."""
    return re.sub(r"(?<=\d)[A-Z](?=~|$)", "", hap)


def read_freqs(path: Path) -> dict[str, float]:
    out = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                hap, f = line.split("\t")[:2]
                out[hap] = float(f)
    return out


def parse_haplotype(hap: str, skip: set[str]) -> dict[str, str]:
    """'A*01:01~C*07:01~B*?' -> {locus: allele}, without unknown alleles and without the loci in `skip`."""
    out = {}
    for part in hap.split("~"):
        locus = part.split("*")[0]
        if "*" in part and not part.endswith("*?") and locus not in skip:
            out[locus] = part
    return out


def family_truth(folder: Path) -> dict[str, tuple[dict, dict]]:
    """sample -> its two haplotypes {locus: allele} from a `hlahaplotools families` output folder."""
    phases = {r["family"]: r for r in read_tsv(folder / "phases.tsv")}
    who = []
    for fam, r in phases.items():
        if r["father"]:
            who.append((r["father"], fam, "A", "B"))
        if r["mother"]:
            who.append((r["mother"], fam, "C", "D"))
    who += [(r["child"], r["family"], r["paternal"], r["maternal"]) for r in read_tsv(folder / "children.tsv")]
    out = {}
    for sample, fam, x, y in who:
        if "?" in (x, y):
            continue
        r = phases[fam]
        bad = set(filter(None, (r["inconsistent_loci"] + "," + r["ambiguous_loci"]).split(",")))

        h1, h2 = parse_haplotype(r[f"hap_{x}"], bad), parse_haplotype(r[f"hap_{y}"], bad)
        if sample not in out or len(set(h1) & set(h2)) > len(set(out[sample][0]) & set(out[sample][1])):
            out[sample] = (h1, h2)
    return out


def pairs(per):
    het = [k for k, (a, b) in enumerate(per) if a != b]
    if not het:
        hap = canon([a for a, _ in per])
        return [(hap, hap)], 0
    out = []
    for flips in itertools.product((False, True), repeat=len(het) - 1):
        flip = dict(zip(het[1:], flips, strict=True))
        out.append((canon([b if flip.get(k) else a for k, (a, b) in enumerate(per)]),
                    canon([a if flip.get(k) else b for k, (a, b) in enumerate(per)])))  # fmt: skip
    return out, len(het)


def cmd_score(args) -> int:
    truths = [family_truth(f) for f in args.families]
    hs_dir = args.haplostats or args.bench / "haplostats"
    summary, people = [], []
    for case in sorted(p for p in (args.bench / "cases").iterdir() if (p / "input.dat").exists()):
        lines = (case / "input.dat").read_text().splitlines()
        loci = lines[0].split("\t")[1::2]
        rows = {c[0]: [(f"{l}*{c[1 + 2 * k]}", f"{l}*{c[2 + 2 * k]}") for k, l in enumerate(loci)]
                for c in (line.split("\t") for line in lines[1:])}  # fmt: skip
        tools = {f.stem.removeprefix("freqs_"): read_freqs(f) for f in case.glob("freqs_*.tsv")}
        if (hs_dir / case.name / "freqs_haplo_stats.tsv").exists():
            tools["haplo_stats"] = read_freqs(hs_dir / case.name / "freqs_haplo_stats.tsv")
        stats = defaultdict(lambda: defaultdict(float))
        n_truth = n_inform = 0
        for sid, per in rows.items():
            truth = next((t[sid] for t in truths if sid in t and all(l in t[sid][0] and l in t[sid][1] for l in loci)),
                         None)  # fmt: skip
            if truth is None:
                continue
            options, nhet = pairs(per)
            # compared without expression suffixes, so a renamed allele (C*04:09N, now C*04:09L) still matches
            true_pair = {bare(canon([truth[0][l] for l in loci])), bare(canon([truth[1][l] for l in loci]))}
            match = next((i for i, (a, b) in enumerate(options) if {bare(a), bare(b)} == true_pair), None)
            n_truth += 1
            if match is None or nhet < 2:
                continue
            n_inform += 1
            for tool, freqs in tools.items():
                scores = [(2 if a != b else 1) * freqs.get(a, 0.0) * freqs.get(b, 0.0) for a, b in options]
                total = sum(scores)
                greater = sum(s > scores[match] for s in scores)
                tied = sum(s == scores[match] for s in scores)
                s = stats[tool]
                s["first"] += (1 / tied) if greater == 0 and scores[match] > 0 else 0
                s["prob"] += scores[match] / total if total > 0 else 0.0
                s["rank"] += 1 + greater + (tied - 1) / 2
                s["zero"] += scores[match] == 0
                people.append([case.name, sid, tool, nhet, len(options), 1 + greater + (tied - 1) / 2,
                               round(scores[match] / total, 4) if total > 0 else 0.0])  # fmt: skip
        for tool, s in sorted(stats.items()):
            summary.append(["~".join(loci), len(rows), n_truth, n_inform, tool,
                            f"{s['first'] / n_inform:.3f}", f"{s['prob'] / n_inform:.3f}",
                            f"{s['rank'] / n_inform:.2f}", int(s["zero"])])  # fmt: skip
    out = args.out or args.bench
    out.mkdir(parents=True, exist_ok=True)
    for name, head, data in (
        ("segregation_check.tsv", ["loci", "people", "with_family_phase", "informative", "tool",
                                   "true_pair_ranked_first", "mean_probability_of_true_pair",
                                   "mean_rank_of_true_pair", "true_pair_given_zero"], summary),
        ("segregation_people.tsv", ["case", "sample", "tool", "heterozygous_loci", "possible_pairs",
                                    "rank_of_true_pair", "probability_of_true_pair"], people),
    ):  # fmt: skip
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(head)
            w.writerows(data)
    for r in summary:
        print("\t".join(map(str, r)))
    print(f"Wrote {out / 'segregation_check.tsv'} and {out / 'segregation_people.tsv'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    f = sub.add_parser("founders", help="write one set of unrelated founders per family as a MAC table")
    f.add_argument("table", type=Path, help="family table (.csv or .tsv)")
    f.add_argument("-o", "--output", type=Path, required=True)
    f.add_argument("--exclude-flags", type=Path, nargs="*", help="flags.tsv files: leave out families flagged as "
                   "parent_excluded or possible_duplicate")  # fmt: skip
    f.set_defaults(func=cmd_founders)
    s = sub.add_parser("score", help="rank each founder's true pair under each tool's frequencies")
    s.add_argument("--bench", type=Path, required=True, help="benchmark.py run --out folder (with cases/)")
    s.add_argument("--families", type=Path, nargs="+", required=True,
                   help="hlahaplotools families output folder(s); the first that phases a founder at every locus "
                   "of a set is used")  # fmt: skip
    s.add_argument("--haplostats", type=Path, help="bench_haplostats.py --out folder (default BENCH/haplostats)")
    s.add_argument("--out", type=Path, help="output folder (default: the --bench folder)")
    s.set_defaults(func=cmd_score)
    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
