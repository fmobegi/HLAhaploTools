"""Plasmode simulation: accuracy at the size, linkage and missing-data pattern of a real cohort.

The haplotype frequencies estimated from the real cohort are taken as the truth. Each
simulated person draws two haplotypes from them, and copies the typed/untyped status of
one real person at each classical locus, so the simulated cohort has the real number of
people, the real haplotype distribution (with its linkage and long rare tail) and the real
pattern of missing loci. DRB3/4/5 count as typed whenever DRB1 is: in the real data a
blank there means an absent gene, which is written as NNNN here. Typings are exact
(two-field names, estimated as typed), so the error measured is that of phase and missing
loci, not of typing ambiguity.

Each replicate is estimated twice with `hlahaplotools estimate`: from everyone, and from
the fully typed people only (all that complete-case tools such as Hapl-o-Mat and PyPop can
use). Error is the total variation distance (TVD) to the truth.

    uv run python scripts/plasmode.py results/full_cohort/classical8/hfs.dat \\
        results/10July2025/haplomap_input_current_names.dat -o bench/plasmode/classical8 --reps 3

Writes summary.json in the output folder.
"""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
import time
from pathlib import Path

CLASSICAL = ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"]
DRB345 = ["DRB3", "DRB4", "DRB5"]


def read_truth(path: Path) -> tuple[list[str], list[list[str]], list[float]]:
    """Loci (in haplotype order), haplotypes as allele lists, frequencies (renormalised)."""
    haps, freqs = [], []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        h, f = line.split("\t")[:2]
        haps.append(h.split("~"))
        freqs.append(float(f))
    loci = [a.split("*")[0] for a in haps[0]]
    total = sum(freqs)
    return loci, haps, [f / total for f in freqs]


def real_patterns(path: Path, loci: list[str], exclude: set[str] = frozenset()) -> list[dict[str, bool]]:
    """Per real person typed at any of `loci` (and not in `exclude`): {locus: typed}."""
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    cols = {l: head.index(l) for l in set(loci) & set(head)}
    out = []
    for line in lines[1:]:
        c = line.split("\t") + [""] * len(head)
        if c[0].strip() in exclude:
            continue
        typed = {l: bool(c[i].strip() or c[i + 1].strip()) for l, i in cols.items()}
        for l in DRB345:
            if l in loci:
                typed[l] = typed.get("DRB1", False)
        if any(typed.get(l, False) for l in loci):
            out.append(typed)
    return out


def write_mac(path: Path, loci: list[str], people: list[list[tuple[str, str] | None]]) -> None:
    with path.open("w") as f:
        f.write("id\t" + "\t".join(l for l in loci for _ in range(2)) + "\n")
        for i, per in enumerate(people):
            cells = []
            for g in per:
                cells += ["", ""] if g is None else [a.split("*", 1)[1] for a in g]
            f.write(f"p{i}\t" + "\t".join(cells) + "\n")


def read_hfs(path: Path) -> dict[str, float]:
    out = {}
    for line in path.read_text().splitlines():
        if line.strip():
            h, f = line.split("\t")[:2]
            out[h] = float(f)
    total = sum(out.values())
    return {h: f / total for h, f in out.items()}


def tvd(p: dict[str, float], q: dict[str, float]) -> float:
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in set(p) | set(q))


COMMON = 1e-3  # "common" haplotypes: true frequency at least 0.1%


def common_error(est: dict[str, float], truth: dict[str, float]) -> dict[str, float]:
    """Error on the common haplotypes, where frequencies matter most in practice.

    tvd_common: half the summed absolute error over haplotypes with true frequency >= COMMON;
    rel_error_common: their mean relative error |est - true| / true.
    """
    common = [k for k, f in truth.items() if f >= COMMON]
    diffs = [abs(est.get(k, 0.0) - truth[k]) for k in common]
    return {
        "tvd_common": round(0.5 * sum(diffs), 5),
        "rel_error_common": round(sum(d / truth[k] for d, k in zip(diffs, common, strict=True)) / len(common), 4),
        "n_common": len(common),
        "truth_mass_common": round(sum(truth[k] for k in common), 4),
    }


def estimate(mac: Path, loci: list[str], out: Path, starts: int) -> dict:
    spec = ",".join(f"{l}:asItIs" for l in loci)
    t = time.perf_counter()
    run = subprocess.run(
        [sys.executable, "-m", "hlahaplotools", "estimate", str(mac), "--loci", spec, "-o", str(out),
         "--starts", str(starts), "--half-typed", "unknown", "--no-drb345", "-q"],
        capture_output=True, text=True,
    )  # fmt: skip
    if run.returncode not in (0, 2):
        raise SystemExit(f"estimate failed for {mac}:\n{run.stderr[-2000:]}")
    info = json.loads((out / "run.json").read_text())
    return {"seconds": round(time.perf_counter() - t, 1), "people": info["samples_used"],
            "incomplete": info["samples_incomplete"], "candidates": info["candidate_haplotypes"],
            "converged": info["converged"]}  # fmt: skip


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("truth", type=Path, help="hfs.dat estimated from the real cohort")
    ap.add_argument("cohort", type=Path, help="the real MAC file (for the missing-locus pattern)")
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--starts", type=int, default=10)
    ap.add_argument("--rescore", action="store_true", help="only recompute the errors of an existing run")
    args = ap.parse_args()

    if args.rescore:
        _, haps, freqs = read_truth(args.truth)
        truth = {"~".join(h): f for h, f in zip(haps, freqs, strict=True)}
        s = json.loads((args.output / "summary.json").read_text())
        for rec in s["reps"]:
            for name in ("everyone", "complete"):
                hfs = args.output / f"rep{rec['rep']}" / name / "hfs.dat"
                if name in rec and hfs.exists():
                    est = read_hfs(hfs)
                    rec[name]["tvd"] = round(tvd(est, truth), 5)
                    rec[name].update(common_error(est, truth))
        (args.output / "summary.json").write_text(json.dumps(s, indent=1) + "\n")
        print(json.dumps([{k: r[k] for k in ("rep", "everyone", "complete") if k in r} for r in s["reps"]], indent=1))
        return 0

    loci, haps, freqs = read_truth(args.truth)
    truth = {"~".join(h): f for h, f in zip(haps, freqs, strict=True)}
    # People the real estimate left out (unreadable typings) are left out here too.
    real_run = args.truth.with_name("run.json")
    left_out = set()
    if real_run.exists():
        info = json.loads(real_run.read_text())
        left_out = {x["id"] for x in info.get("samples_left_out", []) + info.get("unreadable_lines", []) if x.get("id")}
    patterns = real_patterns(args.cohort, loci, left_out)
    print(f"{len(patterns):,} people ({len(left_out)} left out as in the real estimate), {len(haps):,} true haplotypes")
    classical = [l for l in loci if l in CLASSICAL]
    args.output.mkdir(parents=True, exist_ok=True)
    reps = []
    for r in range(args.reps):
        rng = random.Random(args.seed + r)
        people, complete = [], []
        for typed in patterns:
            h1, h2 = rng.choices(haps, weights=freqs, k=2)
            per = [(a, b) if typed.get(l, False) else None for l, a, b in zip(loci, h1, h2, strict=True)]
            people.append(per)
            if all(typed.get(l, False) for l in classical):
                complete.append(per)
        d = args.output / f"rep{r + 1}"
        d.mkdir(exist_ok=True)
        write_mac(d / "everyone.dat", loci, people)
        write_mac(d / "complete.dat", loci, complete)
        rec = {"rep": r + 1, "seed": args.seed + r, "people": len(people), "fully_typed": len(complete)}
        for name in ("everyone", "complete"):
            rec[name] = estimate(d / f"{name}.dat", loci, d / name, args.starts)
            est = read_hfs(d / name / "hfs.dat")
            rec[name]["tvd"] = round(tvd(est, truth), 5)
            rec[name].update(common_error(est, truth))
            print(f"rep {r + 1} {name}: TVD {rec[name]['tvd']}, {rec[name]['seconds']} s", flush=True)
        reps.append(rec)
        (args.output / "summary.json").write_text(json.dumps(
            {"truth": str(args.truth), "cohort": str(args.cohort), "loci": loci, "truth_haplotypes": len(haps),
             "starts": args.starts, "reps": reps}, indent=1) + "\n")  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
