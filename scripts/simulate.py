"""Simulate HLA genotype data from known haplotype frequencies.

Used to (1) compare with Hapl-o-Mat on data with real linkage and (2) check
that estimates with untyped loci stay close to the truth.

    uv run python scripts/simulate.py out/ --samples 5000 --loci A,B,C,DRB1,DQB1 \
        --missing DQB1:0.2,C:0.1 --blank-homozygous 1 --format mac --seed 1

Writes:
    out/population.dat (MAC) or out/population.glc (GLSC)
    out/truth.tsv      true haplotype frequencies (2-field names)

Haplotypes are built from real 2-field allele names in the nomenclature
tables, with Zipf-like frequencies so a few are common and many are rare.
Only the standard library is used.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def two_field_alleles(nomenclature: Path, locus: str) -> list[str]:
    """2-field typings at a locus that expand to at least one expressed allele."""
    out = []
    with (nomenclature / "AllAllelesExpanded.txt").open() as f:
        for line in f:
            key, *alleles = line.rstrip("\n").split("\t")
            if not (key.startswith(locus + "*") and key.count(":") == 1 and key[-1].isdigit()):
                continue
            # Skip names covering only null or unusually expressed alleles (N, L, S, Q).
            if any(a[-1].isdigit() for a in alleles):
                out.append(key)
    return out


def make_haplotypes(rng: random.Random, alleles: dict[str, list[str]], loci: list[str], n_haps: int, n_alleles: int):
    """n_haps distinct haplotypes over a limited allele pool per locus, Zipf weights."""
    pool = {l: rng.sample(alleles[l], min(n_alleles, len(alleles[l]))) for l in loci}
    # Common alleles get picked more often, giving realistic allele skew.
    weights = {l: [1 / (i + 1) for i in range(len(pool[l]))] for l in loci}
    haps: list[tuple[str, ...]] = []
    seen = set()
    while len(haps) < n_haps:
        h = tuple(rng.choices(pool[l], weights[l])[0] for l in loci)
        if h not in seen:
            seen.add(h)
            haps.append(h)
    raw = [1 / (i + 1) ** 1.1 for i in range(n_haps)]
    total = sum(raw)
    return haps, [r / total for r in raw]


# Chromosome 6 order, telomere to centromere (as in hlahaplotools' audit).
GENOMIC_ORDER = [
    "HFE",
    "F",
    "V",
    "G",
    "H",
    "K",
    "A",
    "J",
    "L",
    "E",
    "C",
    "B",
    "MICA",
    "MICB",
    "DRA",
    "DRB9",
    "DRB5",
    "DRB4",
    "DRB3",
    "DRB1",
    "DQA1",
    "DQB1",
    "DQA2",
    "DQB2",
    "DOB",
    "TAP2",
    "TAP1",
    "DMB",
    "DMA",
    "DOA",
    "DPA1",
    "DPB1",
]
ERROR_KINDS = ["other_father", "duplicate_of_mother", "recombinant", "typo"]


def write_families(args, rng: random.Random, loci, haps, freqs, missing: dict[str, float]) -> None:
    """Families: two founders drawn from the haplotype frequencies, children
    inheriting one haplotype from each parent.

    With --errors, some children carry a planted problem:
      other_father         paternal haplotype from an unrelated man (non-paternity / half-sibling)
      duplicate_of_mother  the row is a copy of the mother's typing
      recombinant          paternal haplotype is a crossover of the father's two
      typo                 one allele replaced by another allele of the same locus

    Writes families.csv (FAMILY_ID, Family_Member, LOCUS_1, LOCUS_2, ...),
    truth_founders.tsv (family, F/M, haplotype 1, haplotype 2) and
    truth_errors.tsv (family, member, planted problem).
    """
    lo, _, hi = args.children.partition("-")
    lo, hi = int(lo), int(hi or lo)
    rates = {k: float(v) for k, v in (x.split(":") for x in args.errors.split(",") if x)}
    for k in rates:
        if k not in ERROR_KINDS:
            raise SystemExit(f"unknown error kind {k!r}; use {', '.join(ERROR_KINDS)}")
    order = sorted(range(len(loci)), key=lambda i: GENOMIC_ORDER.index(loci[i]) if loci[i] in GENOMIC_ORDER else 99)
    pool = {l: sorted({h[k] for h in haps}) for k, l in enumerate(loci)}

    def recombine(h1, h2):
        cut = rng.randint(1, len(order) - 1)  # crossover between two loci in chromosome order
        first = set(order[:cut])
        return tuple(h1[k] if k in first else h2[k] for k in range(len(loci)))

    header = ["FAMILY_ID", "Family_Member"] + [f"{l}_{i}" for l in loci for i in (1, 2)]
    table = args.out / "families.csv"
    n_people = 0
    with (
        table.open("w") as out,
        (args.out / "truth_founders.tsv").open("w") as truth,
        (args.out / "truth_errors.tsv").open("w") as errors,
    ):
        out.write(",".join(header) + "\n")
        for fam in range(1, args.samples + 1):
            fid = f"fam{fam:04d}"
            parents = {"F": rng.choices(haps, freqs, k=2), "M": rng.choices(haps, freqs, k=2)}
            members = []
            for role, (h1, h2) in parents.items():
                truth.write(f"{fid}\t{role}\t{'~'.join(h1)}\t{'~'.join(h2)}\n")
                if rng.random() >= args.missing_parent:
                    members.append((role, list(h1), list(h2)))
            for c in range(1, rng.randint(lo, hi) + 1):
                pat, mat = rng.choice(parents["F"]), rng.choice(parents["M"])
                planted = next((k for k in ERROR_KINDS if rng.random() < rates.get(k, 0.0)), None)
                if planted == "other_father":
                    pat = rng.choices(haps, freqs)[0]
                elif planted == "recombinant":
                    # Only keep crossovers that change something (differ from both haplotypes).
                    for _ in range(20):
                        rec = recombine(*rng.sample(parents["F"], 2))
                        if rec not in parents["F"]:
                            pat = rec
                            break
                    else:
                        planted = None
                row = (list(pat), list(mat))
                if planted == "duplicate_of_mother":
                    row = tuple(list(h) for h in parents["M"])
                elif planted == "typo":
                    k = rng.randrange(len(loci))
                    row[0][k] = rng.choice([a for a in pool[loci[k]] if a not in (row[0][k], row[1][k])])
                if planted:
                    errors.write(f"{fid}\t{fid}/C{c}\t{planted}\n")
                members.append((f"C{c}", *row))
            for role, h1, h2 in members:
                cells = []
                for k, locus in enumerate(loci):
                    if rng.random() < missing.get(locus, 0.0):
                        cells += ["", ""]
                    elif h1[k] == h2[k] and rng.random() < args.blank_homozygous:
                        cells += [h1[k], ""]
                    else:
                        cells += [h1[k], h2[k]]
                out.write(",".join([fid, role, *cells]) + "\n")
                n_people += 1
    print(f"wrote {table} ({args.samples} families, {n_people} people), truth_founders.tsv, truth_errors.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--samples", type=int, default=2000)
    ap.add_argument("--loci", default="A,B,C,DRB1,DQB1")
    ap.add_argument("--haplotypes", type=int, default=300, help="distinct true haplotypes")
    ap.add_argument("--alleles", type=int, default=25, help="alleles per locus to draw from")
    ap.add_argument("--missing", default="", help="per-locus untyped rate, e.g. DQB1:0.2,C:0.1")
    ap.add_argument(
        "--half",
        default="",
        help="per-locus rate of one allele dropped at random, e.g. C:0.05 (read with --half-typed unknown)",
    )
    ap.add_argument(
        "--blank-homozygous",
        type=float,
        default=0.0,
        help="chance a homozygous locus is written with the second allele blank, as labs often report it",
    )
    ap.add_argument(
        "--ambiguity",
        type=float,
        default=0.0,
        help="GLSC only: chance a typed allele is reported as a list with 1-2 decoys",
    )
    ap.add_argument(
        "--format",
        choices=["mac", "glsc", "families"],
        default="mac",
        help="families: --samples families of two parents and --children children, as a CSV table",
    )
    ap.add_argument("--children", default="1-4", help="families: children per family, e.g. 2 or 1-4")
    ap.add_argument("--missing-parent", type=float, default=0.0, help="families: chance each parent is not typed")
    ap.add_argument(
        "--errors",
        default="",
        help="families: per-child chance of a planted problem, e.g. "
        "other_father:0.05,duplicate_of_mother:0.02,recombinant:0.05,typo:0.03",
    )
    ap.add_argument("--nomenclature", type=Path, help="default: data/derived/<release from manifest>")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    if args.nomenclature is None:
        import json

        release = json.loads((ROOT / "data" / "manifest.json").read_text())["imgthla"]["release"]
        args.nomenclature = ROOT / "data" / "derived" / release
    rng = random.Random(args.seed)
    loci = args.loci.split(",")
    rate = lambda spec: {k: float(v) for k, v in (x.split(":") for x in spec.split(",") if x)}  # noqa: E731
    missing, half = rate(args.missing), rate(args.half)

    alleles = {l: two_field_alleles(args.nomenclature, l) for l in loci}
    for l, a in alleles.items():
        if not a:
            raise SystemExit(f"no alleles for locus {l} in {args.nomenclature}")
    haps, freqs = make_haplotypes(rng, alleles, loci, args.haplotypes, args.alleles)

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "truth.tsv").open("w") as f:
        for h, p in sorted(zip(haps, freqs, strict=True), key=lambda x: -x[1]):
            f.write("~".join(h) + f"\t{p:.12f}\n")

    if args.format == "families":
        write_families(args, rng, loci, haps, freqs, missing)
        return

    rows, full_rows = [], []
    for i in range(args.samples):
        h1, h2 = rng.choices(haps, freqs, k=2)
        per_locus = []
        for k, locus in enumerate(loci):
            a, b = h1[k], h2[k]
            if rng.random() < missing.get(locus, 0.0):
                per_locus.append(None)
            elif rng.random() < half.get(locus, 0.0):
                per_locus.append((a, None) if rng.random() < 0.5 else (b, None))
            elif a == b and rng.random() < args.blank_homozygous:
                per_locus.append((a, None))
            else:
                per_locus.append((a, b))
        rows.append((f"s{i + 1}", per_locus))
        full_rows.append((f"s{i + 1}", [(h1[k], h2[k]) for k in range(len(loci))]))

    pool = {l: alleles[l] for l in loci}

    def typed(allele: str, locus: str) -> str:
        if rng.random() >= args.ambiguity:
            return allele
        decoys = rng.sample(pool[locus], 2)[: rng.randint(1, 2)]
        items = [allele] + [d for d in decoys if d != allele]
        rng.shuffle(items)
        return "/".join(items)

    def write(path: Path, rows) -> None:
        with path.open("w") as f:
            if args.format == "mac":
                f.write("id\t" + "\t".join(f"{l}\t{l}" for l in loci) + "\n")
                for sid, per_locus in rows:
                    cells = []
                    for g in per_locus:
                        a, b = g if g else (None, None)
                        cells += [a.split("*")[1] if a else "", b.split("*")[1] if b else ""]
                    f.write(sid + "\t" + "\t".join(cells) + "\n")
            else:
                for sid, per_locus in rows:
                    fields = []
                    for locus, g in zip(loci, per_locus, strict=True):
                        if g is None:
                            continue  # untyped loci are simply left out
                        a, b = g
                        fields.append(typed(a, locus) + "+" + (typed(b, locus) if b else ""))
                    f.write(sid + "\t" + "\t".join(fields) + "\n")

    ext = "dat" if args.format == "mac" else "glc"
    path = args.out / f"population.{ext}"
    write(path, rows)
    if missing or half or args.blank_homozygous:
        # The same samples before any typing was removed, for comparison.
        write(args.out / f"population_full.{ext}", full_rows)
        # Only the samples with every locus fully typed (what Hapl-o-Mat would use).
        complete = [r for r in rows if all(g is not None and g[1] is not None for g in r[1])]
        write(args.out / f"population_complete_only.{ext}", complete)
    print(f"wrote {path} ({args.samples} samples) and {args.out / 'truth.tsv'} ({len(haps)} haplotypes)")


if __name__ == "__main__":
    main()
